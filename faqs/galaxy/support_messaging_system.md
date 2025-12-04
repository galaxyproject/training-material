---
title: How can I use the Matrix messaging system for Galaxy Project communication?
area: support
description: When you are directed to use Matrix for a Galaxy working group or for your support question, but you don't know what Matrix or how to access it. Learn what Matrix is, why it is used, and how to use Matrix to connect to us.
box_type: tip
layout: faq
contributors: [hujambo-dunia]
---

### Introduction

The Galaxy community uses [Matrix](https://matrix.org/) – a secure, decentralized, open-source messaging protocol – as its primary chat system. Matrix enables real-time communication across Galaxy contributors, developers, trainers, and users. Common uses include:

* Getting help from the community
* Discussing tool and workflow development
* Organizing training events and materials
* Collaborating on Galaxy infrastructure or scientific questions

Matrix provides an open alternative to proprietary platforms like Slack and Microsoft Teams, with broader flexibility and user privacy in mind.

---

### 🔰 For Newbies: Getting Started with Matrix

If you are new to messaging systems, here's how to get started with Galaxy on Matrix:

<img src="{% link faqs/galaxy/images/support_matrix_element.jpg %}" alt="Galaxy Project Matrix space using the Element client app" width="350" align="right">

#### Step 1: Choose and install a Matrix app

**Install a Matrix client app:** Download and install a Matrix client:

* [**Element**](https://element.io/) – Recommended. Available on [Web](https://app.element.io), [iOS](https://apps.apple.com/app/element/id1083446067), [Android](https://play.google.com/store/apps/details?id=im.vector.app), [Windows/macOS/Linux](https://element.io/download).&#x20;

* [**Cinny**](https://cinny.in/)

* [**Other Matrix clients...**](https://matrix.org/ecosystem/clients/)

All of these apps connect to the same Galaxy channels/rooms. Choose the interface that works best for you.

#### Step 2: Create a Matrix account

* Open the app and sign up. It's easiest to use the default public server at `matrix.org` for your account.
* Choose a username. Your Matrix ID will look like: `@yourname:matrix.org`.

#### Step 3: Join the Galaxy **"Lobby"**

* In your Matrix client, join the [main Galaxy chat room - The Lobby](https://matrix.to/#/#galaxyproject_Lobby:gitter.im).
  * If you're using the Element client, click Explore (the compass icon) and search for "Galaxy" – look for a room named Galaxy or Galaxy Lobby (often with an address like `#galaxyproject:matrix.org`).
  * Alternatively, you can use a direct Matrix link if provided, which will prompt your app to join the room.

* In addition to the [Lobby room](https://matrix.to/#/#galaxyproject_Lobby:gitter.im), join other Matrix rooms via direct links once logged in:

<div align="center">

  <table border="1">
    <thead>
      <tr>
        <th>Galaxy Subjects on Matrix:</th>
        <th>Galaxy Users on Matrix:</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>
          <a href="https://matrix.to/#/#GalaxyProteomics_Lobby:gitter.im" >Proteomics</a> <br>
          <a href="https://matrix.to/#/#galaxyproject_microGalaxy:gitter.im" >Microbiology</a> <br>
          <a href="https://matrix.to/#/#Galaxy-Training-Network_galaxy-single-cell:gitter.im" >Single-Cell & Spatial Omics Users</a> <br>
          <strong><a href="https://matrix.to/#/#galaxyproject:matrix.org" >All Galaxy Matrix rooms, 70+ rooms</a></strong>
        </td>
        <td>
          <a href="https://matrix.to/#/#galaxyproject_tools:gitter.im" >Tool Authors</a> <br>
          <a href="https://matrix.to/#/#galaxyproject_dev:gitter.im" >Developers</a> <br>
          <a href="https://matrix.to/#/#galaxyproject_wg-goat:gitter.im" >Outreach</a> <br>
          <a href="https://matrix.to/#/#galaxyproject_admins:gitter.im" >Admins</a>
        </td>
      </tr>
    </tbody>
  </table>

</div>

**Say hello and ask questions:** Once in the Galaxy Lobby, feel free to introduce yourself or ask your question. This Lobby room is a friendly starting point – community members will welcome you, answer basic questions, and guide you to more specific Galaxy channels/rooms if needed. You will be redirected to the right room if needed.

---

### 💡 For Experienced Users (Slack/Teams/Discord Users)

If you've used tools like Slack, Microsoft Teams, or Discord, Matrix will feel familiar — with a few important differences:

#### ✅ Key Similarities

* **Rooms ≈ Channels**: Rooms are topic-based, like Slack channels.
* **DMs supported**: You can message users privately.
* **Multiple devices**: Stay logged in on phone, tablet, laptop.

#### ❗ Key Differences

* **Federated, decentralized network**: There is no single “Galaxy workspace” to be invited to. Join any Matrix room directly via a public server.
* **Flexible clients**: You can use different Matrix client apps across platforms. They all show the same chats.
* **Room discovery**: Galaxy may provide a **Matrix Space** to group related rooms to organize communities (e.g. training, user help, development, working groups/wg, infrastructure).
* **Bridged rooms**: Some rooms are connected to Gitter or other services, but behave normally in Matrix. For example:

  ```
  #galaxyproject_Lobby:gitter.im
  ```

  This is a Matrix room that also syncs with users on Gitter. You can treat them as normal Matrix rooms; messages are synced across the platforms.

---

### 🔐 Security and Privacy

Matrix supports **end-to-end encryption (E2EE)** in private conversations and invite-only rooms.

However, **public Galaxy rooms are not encrypted**. This is intentional so people can easily join and search history. Therefore:

* 🔍 **Assume public visibility**: Don't share passwords, private research data, or anything sensitive.
* 🙈 **Use nicknames if preferred**: Your Matrix ID is visible in public rooms, but you don't have to use your real name.
* 🔐 **Private DMs are encrypted by default**.
* ⚠️ **Be cautious**: Even with encryption, nothing is foolproof. Use common sense.

Matrix is designed for open collaboration. When in doubt, treat public rooms like an open forum.

---

### 📍 TL;DR Quickstart

* ✅ Install [Element](https://element.io) or another Matrix app
* ✅ Create an account (Matrix ID)
* ✅ Join `#galaxyproject:matrix.org` (Galaxy Lobby)
* ✅ Ask questions or join other Galaxy rooms from there
* ❌ Don't share private info in public rooms

For more Matrix help: [https://matrix.org/docs](https://matrix.org/docs)

See you in the Lobby! 🎉
