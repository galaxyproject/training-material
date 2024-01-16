---
title: Blank page or no CSS/JavaScript
area: deployment
box_type: tip
layout: faq
contributors: [natefoo]
---

This generally means that serving of static content is broken:

- Check browser console for 404 errors.
- Check proxy error log for permission errors.
- Verify that your proxy static configuration is correct.
- If you have recently upgraded Galaxy or changed the GUI in some way, you will need to rebuild the client
