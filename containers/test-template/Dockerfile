ARG PYTHON_VERSION=3.11
FROM python:${PYTHON_VERSION}-slim

# Install and uninstall snakebids to cache it and it's dependences
RUN apt-get update && apt-get install -y gcc git && \
    rm -rf /var/lib/apt/lists/* && \
    git config --global --add safe.directory /src && \
    python -m pip install pipx && \
    pipx install poetry && \
    pipx install hatch && \
    pipx install pdm && \
    mkdir prebuild && \
    cd prebuild && \
    pip wheel snakebids && \
    cd .. && \
    rm -rf prebuild

COPY ./test-template.sh /run/test-template.sh
ENV PATH="/root/.local/bin:$PATH"

WORKDIR /work
ENTRYPOINT [ "/run/test-template.sh" ]
