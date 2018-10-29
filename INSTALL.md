# Installation Requirements

<!-- TOC depthFrom:2 -->

- [Debian](#debian)
- [RHEL/CentOS](#rhelcentos)

<!-- /TOC -->

## Debian

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    autoconf \
    pkg-config \
    zlib1g-dev \
    check
```

## RHEL/CentOS

```bash
sudo yum groupinstall -y "Development Tools"
sudo yum install -y \
    bzip2 \
    zlib-devel \
    check-devel
```
