---
title: "Module 2_2: Setup Data for ‘Introduction to the Command Line for Genomics’"
layout: page
permalink: /tutorials/cli-genomics-setup/
date: 2025-05-21
description: "Download and unpack the Data Carpentry shell_data files on your Vertex AI Workbench (or any Linux terminal)."
nav_order: 3
---

> These instructions prepare the **shell_data** directory used in the Software Carpentry *Introduction to the Command Line for Genomics* lesson. Follow them in any Linux terminal—including your Vertex AI Workbench instance.

---

## 1 · Start in your home directory
Always keep workshop data in a tidy location (e.g., `$HOME` or a `data/` sub-folder).

```bash
cd ~    # jump to your home directory
```

---

## 2 · Download the dataset archive
Use **`wget`** to fetch the compressed file from Figshare (≈ 150 MB).

```bash
wget https://figshare.com/ndownloader/articles/7726454/versions/3
```

*Result:* a file called **`3`** appears in your current directory—it’s actually a ZIP archive.

---

## 3 · Rename, extract, and unpack

1. **Rename** the file for clarity:
   ```bash
   mv 3 cli.zip
   ```

2. **Extract** the `shell_data.tar.gz` tarball from the ZIP (no need to unzip everything else):
   ```bash
   unzip cli.zip shell_data.tar.gz
   ```

3. **Unpack** the tarball to create the **shell_data/** directory:
   ```bash
   tar xvf shell_data.tar.gz
   ```

> Commands explained  
> • `mv` – move/rename files  
> • `unzip` – extract from ZIP archives  
> • `tar xvf` – e**x**tract (**v**erbose) the specified **f**ile

---

## 4 · Verify your data
List the contents of the new directory; the trailing `/` tells `ls` you’re looking at a folder.

```bash
ls -F shell_data/
```

You should see sub-folders like `sra_metadata/` and `untrimmed_fastq/`. If so, you’re ready for the lesson!

---

### Next steps
Open the first lesson exercise at https://datacarpentry.github.io/shell-genomics/01-introduction.html
