---
title: "Tutorial: Create Google Cloud Vertex AI Workbench"
layout: page
permalink: /tutorials/vertex-ai-workbench/
date: 2025-05-21
description: "Create your first Vertex AI Workbench instance using **only** the Google Cloud Console."
nav_order: 2
---

> **Why Vertex AI Workbench?** It provides a managed JupyterLab environment that lives close to your Google Cloud data and scales from modest CPU notebooks to GPU/TPU powerhouses—perfect for genomics and machine‑learning workflows.

This quick‑start shows how to spin up a **user‑managed notebook** instance *entirely through the Google Cloud Console*.

## Prerequisites

1. A Google Cloud project **with billing enabled**.
2. **Owner** or **Editor** permissions in that project.

---

## 1 · Enable required APIs (Console)

1. In the left sidebar, select **APIs & Services → Library**.
2. Search for and enable the following:
   - **Vertex AI API**
   - **Notebooks API**
   - **Compute Engine API**

> You’ll see an **Enable** button on each API page—click it, wait a few seconds, then move on to the next API.

---

## 2 · Note a region

Vertex AI Workbench instances are **regional**. If you’re based in Alaska, `us‑west1 (Oregon)` or `us‑central1 (Iowa)` usually give the best latency. Jot down the region you’ll use; you’ll pick it again in the next step.

---

## 3 · Create the Workbench instance

1. Navigate to **Vertex AI → Workbench**.
2. Make sure you're in the **Instances** view.
3. Click **Create New**.
4. Fill in the form:

| Field            | Recommended value                                  |
| ---------------- | -------------------------------------------------- |
| **Name**         | `vertex-genomics-demo-initials`                    |
| **Region**       | The region you chose above (e.g., **us-central1**) |
| **Machine type** | **e2-standard-4** (4 vCPU / 16 GB RAM)             |

5. Leave **Permissions** at the default. (Workbench automatically creates/uses a service account with the *Notebooks Service Agent* role.)
6. Click **Create**. The instance status turns **PROVISIONING**, then **RUNNING**—usually within two minutes.

---

## 4 · Open JupyterLab

1. Once status shows **RUNNING**, click the instance name.
2. Click **Open JupyterLab**. A new browser tab opens your notebook environment—no SSH keys or port forwarding needed.

---

## 5 · Verify the environment with “Hello, world!” examples

### Python notebook

1. In the JupyterLab **Launcher** (the tab that opens by default), click the **Python 3** notebook icon to start a new notebook.  
2. In the first cell, type:

```python
print("Hello, world!")
```

3. Press **Shift + Enter**. You should see `Hello, world!` printed below the cell—confirming the Python kernel is active.

### Bash terminal

1. Still in the **Launcher**, click the **Terminal** icon to open a shell tab.  
2. At the prompt, type:

```bash
echo "Hello, world from Bash!"
```

3. Hit **Enter**. The terminal should echo the same text back, proving the VM’s shell environment is working.

> Feel free to explore additional kernels (R, Julia) or open more terminals from the Launcher.

---

## 6 · Cost‑saving tips

- **Stop** the instance when idle: return to the **Instances** page (Vertex AI → Workbench), tick the checkbox next to your notebook VM, and click **Stop** in the toolbar.
