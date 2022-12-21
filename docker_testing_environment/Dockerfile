FROM alpine:latest AS builder

RUN apk add curl

RUN \
  curl \
  --location https://github.com/grimme-lab/xtb/releases/download/v6.5.1/xtb-6.5.1-linux-x86_64.tar.xz \
  --output xtb.tar.xz \
  && tar xf xtb.tar.xz

FROM python:3.10.8-slim

COPY --from=builder /xtb-6.5.1 /opt/xtb
ENV PATH="/opt/xtb/bin:$PATH"

# requirements.txt is generated with pip-compile.
COPY requirements.txt /deps/requirements.txt

RUN python -m venv /deps/venv && . /deps/venv/bin/activate && pip install -r /deps/requirements.txt
