---
title: Transfer entire histories from one Galaxy server to another
area: histories 
box_type: tip
layout: faq
contributors: [jennaj, AnomalyCodes]
---

**Transfer a Single Dataset**

At the **sender** Galaxy server, [set the history to a shared state]({% link faqs/galaxy/histories_sharing.md %}), then directly capture the {% icon galaxy-link %} link for a dataset and paste the URL into the **Upload** tool at the **receiver** Galaxy server. 

**Transfer an Entire History**

[Have an account]({% link faqs/galaxy/account_create.md %}) at two different Galaxy servers, and be logged into both.

At the **sender** Galaxy server

1. Navigate to the history you want to transfer, and [set the history to a shared state]({% link faqs/galaxy/histories_sharing.md %}).
2. Click into the **History Options** menu in the history panel.
3. Select from the menu {% icon galaxy-history-archive %} **Export History to File**.
4. Choose the option for **How do you want to export this History?** as **to direct download**.
5. Click on **Generate direct download**.
6. Allow the archive generation process to complete. \*
7. Copy the {% icon galaxy-link %} link for your new archive.

At the **receiver** Galaxy server

8. Confirm that you are logged into your account.
9. Click on **Data** in the top menu, and choose **Histories** to reach your **Saved Histories**.
10. Click on **Import history** in the grey button on the top right.
11. Paste in your link's URL from step 7.
12. Click on **Import History**.
13. Allow the archive import process to complete. \*
14. The transfered history will be uncompressed and added to your **Saved Histories**.


\* For steps 6 and 13: It is Ok to navigate away for other tasks during processing. If enabled, Galaxy will send you [status notifications]({% link faqs/galaxy/account_update_preference.md %}).


{% icon tip %} If the history to transfer is large, you may [copy just your important datasets into a new history]({% link faqs/galaxy/histories_copy_dataset.md %}), and create the archive from that new smaller history. Clearing away deleted and purged datasets will make *all* histories smaller and faster to archive and transfer!

