---
title: How do I load data using FTP?
area: Loading data
layout: faq
box_type: tip
contributors: [jennaj, bernandez]
---

1. Load the data using line command FTP or a client. Note that the FTP server name is specific to the Galaxy you are working on. This is by default the URL of the server.
    * For the public Galaxy Main instance at [http://usegalaxy.org](http://usegalaxy.org) the FTP server name to use is **usegalaxy.org**.
    * For a default local (with FTP enabled, see next) the FTP server name to use is **localhost:8080**. If the server URL was modified, use that custom URL.
    * **`FTPS` was enabled for all transfers to [http://usegalaxy.org](http://usegalaxy.org) on July 19, 2017**. If you are having trouble connecting the first      time after this date, verifying the server certificate is required when using an FTP client.
    * Working at [Galaxy EU](https://usegalaxy.eu)? Read the server-specific FTP help [here](https://galaxyproject.eu/ftp/):
    * More help for FTP is at Galaxy [Help](https://help.galaxyproject.org). Search with the keyword "ftp". Example post: https://help.galaxyproject.org/t/ftp-help-guides-tutorials-and-troubleshooting/3449.
2. If on another server, the FTP server name will appear in the **Upload** tool pop-up window. When using a local Galaxy server, be certain to configure your instance for FTP first.
4. Use your email and password for the same instance as your credentials to log in and save the data to your account.
5. Once the data is loaded (confirm through FTP client), use the **Upload** tool to load the data into a History.
   
