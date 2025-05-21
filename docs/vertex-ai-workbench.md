---
title: "Tutorial: Create Google Cloud Vertex AI Workbench"
layout: page
permalink: /tutorials/vertex-ai-workbench/
date: 2025-05-21
description: "Create your first Vertex AI Workbench instance using **only** the Google Cloud Console—no gcloud CLI required."
nav_order: 2
---

> **Why Vertex AI Workbench?** It provides a managed JupyterLab environment that lives close to your Google Cloud data and scales from modest CPU notebooks to GPU/TPU powerhouses—perfect for genomics and machine‑learning workflows.

This quick‑start shows how to spin up a **user‑managed notebook** instance *entirely through the Google Cloud Console*.

## Prerequisites

1. A Google Cloud project **with billing enabled**.
2. **Owner** or **Editor** permissions in that project.
3. A modern web browser.

---

## 1 · Enable required APIs (Console)

1. In the left sidebar, select **APIs & Services → Library**.
2. Search for and enable the following:
   - **Vertex AI API**
   - **Notebooks API**
   - **Compute Engine API**

> You’ll see an **Enable** button on each API page—click it, wait a few seconds, then move on to the next API.

---

## 2 · Choose a region

Vertex AI Workbench instances are **regional**. If you’re based in Alaska, `us‑west1 (Oregon)` or `us‑central1 (Iowa)` usually give the best latency. Jot down the region you’ll use; you’ll pick it again in the next step.

---

## 3 · Create the Workbench instance

1. Navigate to **Vertex AI → Workbench → User‑managed notebooks**.
2. Click **New → Customize**.
3. Fill in the form:

| Field            | Recommended value                               |
| ---------------- | ----------------------------------------------- |
| **Name**         | `vertex‑genomics‑demo`                          |
| **Region**       | The region you chose above (e.g., **us‑west1**) |
| **Machine type** | **n1‑standard‑4** (4 vCPU / 15 GB RAM)          |
| **GPUs**         | *None* (add 1× T4 if you need CUDA)             |
| **Boot disk**    | 100 GB, Debian 11                               |

4. Leave **Permissions** at the default. (Workbench automatically creates/uses a service account with the *Notebooks Service Agent* role.)
5. Click **Create**. The instance status turns **PROVISIONING**, then **RUNNING**—usually within two minutes.

---

## 4 · Open JupyterLab

1. Once status shows **RUNNING**, click the instance name.
2. Click **Open JupyterLab**. A new browser tab opens your notebook environment—no SSH keys or port forwarding needed.

---

## 5 · First‑run housekeeping *(optional)*

Inside JupyterLab, open a **Terminal** and consider updating packages:

```bash
sudo apt-get update && sudo apt-get -y upgrade
pip install --upgrade pip
conda update -y conda
```

*(These commands run inside the notebook VM; they don’t require the gcloud CLI.)*

---

## 6 · Cost‑saving tips

- **Stop** the instance when idle: three‑dot menu → **Stop**.
- Set an **idle shutdown** timer: **Settings → Advanced → Idle timeout**.
- Delete unused disks and snapshots in **Compute Engine → Disks**.

---

## Next steps

| Task                                | Where to find details                  |
| ----------------------------------- | -------------------------------------- |
| Attach a GPU later                  | **Vertex AI → Workbench → Edit**       |
| Mount data from Cloud Storage       | **JupyterLab → Extensions → GCS Fuse** |
| Train at scale with Vertex Training | **Vertex AI → Training**               |

---

### Reference

Full Google guide (console‑only version): <https://cloud.google.com/vertex-ai/docs/workbench/instances/create?hl=en>  
Document last tested **May 2025**.
