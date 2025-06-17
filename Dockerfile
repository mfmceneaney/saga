# Start from ROOT base image
FROM rootproject/root:latest

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install build tools, Doxygen, venv support
RUN apt-get update && \
    apt-get install -y \
        cmake \
        g++ \
        git \
        doxygen \
        python3-venv \
        python3-pip \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create Python virtual environment with system Python
RUN python3 -m venv $VIRTUAL_ENV

# Set working directory
WORKDIR /usr/src/app
COPY . .

# Install pip packages for docs
COPY docs/requirements.txt .
RUN pip install --upgrade pip && \
    if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

# Build project and docs with CMake
RUN cmake -S . -B build -DBUILD_DOXYGEN=TRUE && cmake --build build

# Install python modules
RUN if [ -f pyproject.toml ]; then pip install -e .; fi

CMD ["/bin/bash"]
