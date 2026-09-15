## Building from Source

GNU Autotools, GCC (or Clang), and Make are required to build `ceif`. 

### Dependencies & Prerequisites

`ceif` uses `json-c` (or `libfastjson`) to save and load forest models in JSON format by default. If JSON development headers are missing at build time, `./configure` will issue a warning, and `ceif` will compile with a fallback stub that saves models in CSV format.

To build `ceif` with full JSON support, install the development packages for your distribution:

#### Rocky Linux 9 / RHEL 9 / AlmaLinux 9 / CentOS Stream 9
On RHEL-based enterprise distributions, `json-c-devel` resides in the **CRB** (CodeReady Linux Builder) repository, which is disabled by default:
```bash
sudo dnf --enablerepo=crb install -y json-c-devel gcc make autoconf automake
```

#### Fedora
```bash
sudo dnf install -y json-c-devel gcc make autoconf automake
```

#### Debian / Ubuntu / Linux Mint
```bash
sudo apt-get update
sudo apt-get install -y libjson-c-dev gcc make autoconf automake
```

#### Alpine Linux
```bash
apk add json-c-dev gcc make autoconf automake musl-dev
```

#### Arch Linux / Manjaro
```bash
sudo pacman -S json-c base-devel
```

#### macOS (Homebrew)
```bash
brew install json-c autoconf automake
```

---

### Compiling `ceif`

Clone from GitHub and build:

```bash
cd ceif
autoreconf -is
./configure
make
```

Verify that `./configure` detected the JSON library:
```text
checking for library containing json_object_new_object... -ljson-c
checking for json-c/json.h... yes
```

### Installation & Verification

Test the compiled binary:
```bash
./src/ceif -V
```

Optionally install system-wide (defaults to `/usr/local/bin`):
```bash
sudo make install
```
