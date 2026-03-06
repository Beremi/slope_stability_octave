FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV REQUIRE_QT_TOOLKIT=1

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash-completion \
        bison \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        flex \
        gawk \
        gfortran \
        git \
        gnuplot-nox \
        libarpack2-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libegl-dev \
        libfftw3-dev \
        libfontconfig-dev \
        libgl-dev \
        libgl2ps-dev \
        libglu1-mesa-dev \
        libgraphicsmagick++1-dev \
        libhdf5-dev \
        liblzma-dev \
        libncurses-dev \
        libpcre2-dev \
        libportaudio2 \
        libqscintilla2-qt6-dev \
        qt6-5compat-dev \
        libqt6opengl6-dev \
        libqt6svg6-dev \
        libreadline-dev \
        libsndfile1-dev \
        libz-dev \
        m4 \
        mesa-common-dev \
        ninja-build \
        perl \
        pkg-config \
        python3 \
        python3-pip \
        python3-venv \
        qt6-base-dev \
        qt6-tools-dev \
        qt6-tools-dev-tools \
        sed \
        tar \
        wget \
        xauth \
        xvfb \
        xz-utils \
    && rm -rf /var/lib/apt/lists/*

COPY .devcontainer/auto-activate.sh /usr/local/share/slope-stability/auto-activate.sh

RUN chmod 0644 /usr/local/share/slope-stability/auto-activate.sh \
    && printf '\n# Auto-activate local slope_stability toolchain when present\nsource /usr/local/share/slope-stability/auto-activate.sh\n' >> /etc/bash.bashrc \
    && mkdir -p /workspaces

WORKDIR /workspaces
USER vscode
