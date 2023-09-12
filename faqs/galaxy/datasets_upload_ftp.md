---
title: Upload many files (>10) via FTP
area: data upload
box_type: tip
layout: faq
contributors: [bebatut]
---

<div class="alert alert-success trim-p" role="alert">
UPDATE: FTP upload is no longer supported on UseGalaxy.org. For more information and to learn about alternative options, please see https://help.galaxyproject.org/t/the-usegalaxy-org-ftp-service-will-be-decommissioned-on-august-12-2022/8318.
</div>

1. Make sure to have an FTP client installed

    There are many options. We can recommend [FileZilla](https://filezilla-project.org/), a free FTP client that is available on Windows, MacOS, and Linux.

2. Establish FTP connection to the Galaxy server
    1. Provide the Galaxy server's FTP server name (e.g. `usegalaxy.org`, `ftp.usegalaxy.eu`)
    2. Provide the **username** (usually the email address) and the **password** on the Galaxy server
    3. Connect

3. Add the files to the FTP server by dragging/dropping them or right clicking on them and uploading them

    The FTP transfer will start. We need to wait until they are done.

4. Open the Upload menu on the Galaxy server
5. Click on **Choose FTP file** on the bottom
6. Select files to import into the history
7. Click on **Start**
