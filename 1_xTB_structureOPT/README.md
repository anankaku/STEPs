# Installing xTB via Conda (Linux)

This guide describes how to install **xTB** in a clean Linux environment using **Conda**.
The procedure does **not** require root privileges and is suitable for local machines, HPC clusters, and cloud environments.

---

## System Requirements

* Operating System: Linux (x86_64)
* User privileges: standard user (no `sudo` required)
* Internet access (for initial installation)

---

## Why Use Conda?

* No root access required
* Automatic dependency resolution (OpenMP, C++ runtime, BLAS)
* Reproducible environments
* Compatible with STEPs / Amber / OpenMM workflows

---

## 1. Install Miniconda

Download the Miniconda installer:

```bash
cd ~
```
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Run the installer:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Recommended options during installation:

* Accept the license (`yes`)
* Use the default installation path (`~/miniconda3`)
* Allow Conda to initialize your shell (`yes`)

Reload your shell environment:

```bash
source ~/.bashrc
```

Verify Conda installation:

```bash
conda --version
```

---

## 2. Create a Dedicated xTB Environment

Using a dedicated Conda environment is strongly recommended.

```bash
conda create -n xtb python=3.10
conda activate xtb
```

---

## 3. Configure Conda Channels

Enable `conda-forge` and enforce strict channel priority:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

---

## 4. Install xTB

Install xTB from `conda-forge`:

```bash
conda install xtb
```

This will automatically install all required runtime dependencies.

---

## 5. Verify Installation

### Check version

```bash
xtb --version
```

A successful installation will print the installed xTB version.

---

### Minimal test calculation

Create a simple test system:

```bash
cat > test.xyz << EOF
2
H2 molecule
H  0.0  0.0  0.0
H  0.0  0.0  0.74
EOF
```

Run xTB [*(xTB Method Selection)*](METHOD.md):

```bash
xtb test.xyz --gfn 2
```

If the calculation finishes normally and reports an energy, the installation is working correctly.

---

## Notes

* xTB can be used **standalone** without ORCA.
* ORCA is only required if xTB is invoked as part of an ORCA-based workflow (e.g., mdgx integration).
* On HPC systems, this Conda-based installation does not require administrator intervention.
