FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

ARG VERSION
RUN pip install --upgrade pip && \
    pip install "sniffcell==${VERSION}"

ENTRYPOINT ["sniffcell"]
CMD ["-h"]
