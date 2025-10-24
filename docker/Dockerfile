# Start from ROOT base image #NOTE: If you are using podman to emulate docker then you must use the below ROOT specification form.
FROM docker.io/rootproject/root:6.36.00-ubuntu25.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install build tools, Doxygen, venv support and LaTeX deps for building docs
RUN apt-get update && \
    apt-get install -y \
        bash \
        cmake \
        g++ \
        git \
        doxygen \
        python3-venv \
        python3-pip \
        texlive-latex-recommended \
        texlive-latex-extra \
        texlive-fonts-recommended \
        texlive-xetex \
        latexmk \
        dvipng \
        cm-super \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create Python virtual environment with system Python
RUN python3 -m venv $VIRTUAL_ENV

# Set working directory
WORKDIR /usr/src/saga
COPY . .

# Install pip packages for docs
COPY docs/requirements.txt .
RUN pip install --upgrade pip && \
    if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

# Build project and docs with CMake
RUN cmake -S . -B build -DBUILD_DOXYGEN=TRUE && \
    cmake --build build && \
    cmake --install build --prefix /usr/bin

# Install python modules
RUN if [ -f pyproject.toml ]; then pip install -e .; fi

CMD ["/bin/bash","-c","cd /usr/src/saga && source /usr/src/saga/bin/env.sh && exec bash"]
