classdef Profiler < handle
    %--------------------------------------------------------------------------
    % Profiler  Lightweight wall-clock profiler for class-level instrumentation.
    %
    % A single Profiler handle is shared (by reference) across multiple
    % objects (e.g. CONSTITUTIVE and DFGMRES).  Each instrumented method
    % reports its elapsed time via add_time(); the profiler accumulates
    % totals and call counts automatically.
    %
    % Named entries that contain a '/' are treated as *sub-entries* of their
    % parent (everything before the '/').  print_summary() shows top-level
    % entries sorted by time, then sub-profiles grouped by parent, then
    % counters.
    %
    % Usage:
    %   prof = PROFILING.Profiler();          % create once
    %   my_object.profiler = prof;            % attach (handle semantics)
    %   ...run computation...
    %   prof.print_summary();                 % formatted report
    %   s = prof.to_struct();                 % export for post-processing
    %
    %--------------------------------------------------------------------------

    properties
        data       % containers.Map  name -> struct(total_time, n_calls)
        counters   % containers.Map  name -> double
        enabled    % logical  (false => all calls are no-ops)
    end

    methods

        function obj = Profiler(enabled)
            %--------------------------------------------------------------
            % Profiler  Constructor.
            %   obj = PROFILING.Profiler()          % enabled by default
            %   obj = PROFILING.Profiler(false)     % disabled (zero cost)
            %--------------------------------------------------------------
            if nargin < 1, enabled = true; end
            obj.enabled = enabled;
            obj.data     = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.counters = containers.Map('KeyType', 'char', 'ValueType', 'double');
        end

        function add_time(obj, name, elapsed)
            %--------------------------------------------------------------
            % add_time  Record a timing measurement.
            %   obj.add_time('CONSTITUTIVE.reduction', elapsed)
            %--------------------------------------------------------------
            if ~obj.enabled, return; end
            if obj.data.isKey(name)
                entry = obj.data(name);
                entry.total_time = entry.total_time + elapsed;
                entry.n_calls    = entry.n_calls + 1;
                obj.data(name) = entry;
            else
                obj.data(name) = struct('total_time', elapsed, 'n_calls', 1);
            end
        end

        function add_count(obj, name, n)
            %--------------------------------------------------------------
            % add_count  Increment a named counter (e.g. GMRES iterations).
            %   obj.add_count('DFGMRES.solve.n_gmres_iters', iters)
            %--------------------------------------------------------------
            if ~obj.enabled, return; end
            if nargin < 3, n = 1; end
            if obj.counters.isKey(name)
                obj.counters(name) = obj.counters(name) + n;
            else
                obj.counters(name) = n;
            end
        end

        function print_summary(obj, prefix_filter)
            %--------------------------------------------------------------
            % print_summary  Pretty-print the profiler report.
            %
            %   obj.print_summary()            % full report
            %   obj.print_summary('DFGMRES')   % only DFGMRES entries
            %--------------------------------------------------------------
            if nargin < 2, prefix_filter = ''; end

            all_keys = obj.data.keys();
            if isempty(all_keys)
                fprintf('Profiler: no data recorded.\n');
                return;
            end

            % Optional prefix filter
            if ~isempty(prefix_filter)
                mask = cellfun(@(k) strncmp(k, prefix_filter, length(prefix_filter)), all_keys);
                all_keys = all_keys(mask);
            end
            if isempty(all_keys)
                fprintf('Profiler: no entries matching "%s".\n', prefix_filter);
                return;
            end

            % Separate top-level (no '/') from sub-entries (contain '/')
            is_sub = cellfun(@(k) any(k == '/'), all_keys);
            top_keys = all_keys(~is_sub);
            sub_keys = all_keys(is_sub);

            % ---- Top-level table, sorted by time descending ----
            if ~isempty(top_keys)
                top_times = cellfun(@(k) obj.data(k).total_time, top_keys);
                top_calls = cellfun(@(k) obj.data(k).n_calls,    top_keys);
                [top_times, idx] = sort(top_times, 'descend');
                top_keys  = top_keys(idx);
                top_calls = top_calls(idx);
                total_time = sum(top_times);

                fprintf('\n============================================================\n');
                fprintf('  Profiler Summary\n');
                fprintf('============================================================\n');
                fprintf('  %-9s %-6s %6s  %s\n', 'Time', '%', 'Calls', 'Operation');
                fprintf('  %-9s %-6s %6s  %s\n', '---------', '------', '------', '------------------------------');
                for k = 1:numel(top_keys)
                    pct = 100 * top_times(k) / max(total_time, 1e-12);
                    fprintf('  %8.2fs %5.1f%% %6d  %s\n', top_times(k), pct, top_calls(k), top_keys{k});
                end
                fprintf('  %-9s %-6s %6s  %s\n', '---------', '------', '------', '------------------------------');
                fprintf('  %8.2fs                %s\n', total_time, 'TOTAL');
                fprintf('============================================================\n');
            end

            % ---- Sub-profiles, grouped by parent ----
            if ~isempty(sub_keys)
                % Find unique parents
                parents = cell(size(sub_keys));
                for k = 1:numel(sub_keys)
                    slash_pos = find(sub_keys{k} == '/', 1);
                    parents{k} = sub_keys{k}(1:slash_pos - 1);
                end
                unique_parents = unique(parents);

                for p = 1:numel(unique_parents)
                    parent = unique_parents{p};
                    pmask  = strcmp(parents, parent);
                    pkeys  = sub_keys(pmask);
                    ptimes = cellfun(@(k) obj.data(k).total_time, pkeys);
                    pcalls = cellfun(@(k) obj.data(k).n_calls,    pkeys);
                    [ptimes, pidx] = sort(ptimes, 'descend');
                    pkeys  = pkeys(pidx);
                    pcalls = pcalls(pidx);

                    % Use parent total for %
                    if obj.data.isKey(parent)
                        parent_total = obj.data(parent).total_time;
                        parent_calls = obj.data(parent).n_calls;
                    else
                        parent_total = sum(ptimes);
                        parent_calls = 0;
                    end

                    fprintf('\n  Sub-profile: %s', parent);
                    if parent_calls > 0
                        fprintf('  (%.2fs total, %d calls)', parent_total, parent_calls);
                    end
                    fprintf('\n');
                    fprintf('  %-9s %-6s %6s  %s\n', '---------', '------', '------', '------------------------------');
                    for k = 1:numel(pkeys)
                        child_name = pkeys{k}(length(parent) + 2 : end);
                        pct = 100 * ptimes(k) / max(parent_total, 1e-12);
                        fprintf('  %8.2fs %5.1f%% %6d  %s\n', ptimes(k), pct, pcalls(k), child_name);
                    end
                end
            end

            % ---- Counters ----
            ckeys = obj.counters.keys();
            if ~isempty(prefix_filter)
                cmask = cellfun(@(k) strncmp(k, prefix_filter, length(prefix_filter)), ckeys);
                ckeys = ckeys(cmask);
            end
            if ~isempty(ckeys)
                % Sort alphabetically
                ckeys = sort(ckeys);
                fprintf('\n  Counters:\n');
                for k = 1:numel(ckeys)
                    fprintf('    %-45s = %d\n', ckeys{k}, obj.counters(ckeys{k}));
                end
            end

            fprintf('============================================================\n');
        end

        function s = to_struct(obj)
            %--------------------------------------------------------------
            % to_struct  Export profiler data as a plain struct.
            %--------------------------------------------------------------
            s = struct();
            s.times = struct();
            keys = obj.data.keys();
            for k = 1:numel(keys)
                safe = strrep(strrep(keys{k}, '.', '__'), '/', '___');
                entry = obj.data(keys{k});
                s.times.(safe) = struct('total_time', entry.total_time, ...
                    'n_calls',    entry.n_calls);
            end
            s.counters = struct();
            ckeys = obj.counters.keys();
            for k = 1:numel(ckeys)
                safe = strrep(strrep(ckeys{k}, '.', '__'), '/', '___');
                s.counters.(safe) = obj.counters(ckeys{k});
            end
        end

        function reset(obj)
            %--------------------------------------------------------------
            % reset  Clear all recorded data.
            %--------------------------------------------------------------
            obj.data     = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.counters = containers.Map('KeyType', 'char', 'ValueType', 'double');
        end

    end
end
