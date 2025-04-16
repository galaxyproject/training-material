---
title: Why is my tool erroring as 'Above error raised while reading key '/layers' of type <class 'h5py._hl.group.Group'> from /.'
area: Analysis
box_type: tip
layout: faq
contributors: [nomadscientist]
---

Are you getting the following error, or similar?
```
Traceback (most recent call last):
File "/usr/local/lib/python3.9/site-packages/anndata/_io/utils.py", line 177, in func_wrapper
return func(elem, *args, **kwargs)
File "/usr/local/lib/python3.9/site-packages/anndata/_io/h5ad.py", line 527, in read_group
EncodingVersions[encoding_type].check(
File "/usr/local/lib/python3.9/enum.py", line 432, in __getitem__
return cls._member_map_[name]
KeyError: 'dict'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
File "/usr/local/bin/scanpy-cli", line 10, in <module>
sys.exit(cli())
File "/usr/local/lib/python3.9/site-packages/click/core.py", line 829, in __call__
return self.main(*args, **kwargs)
File "/usr/local/lib/python3.9/site-packages/click/core.py", line 782, in main
rv = self.invoke(ctx)
File "/usr/local/lib/python3.9/site-packages/click/core.py", line 1259, in invoke
return _process_result(sub_ctx.command.invoke(sub_ctx))
File "/usr/local/lib/python3.9/site-packages/click/core.py", line 1259, in invoke
return _process_result(sub_ctx.command.invoke(sub_ctx))
File "/usr/local/lib/python3.9/site-packages/click/core.py", line 1066, in invoke
return ctx.invoke(self.callback, **ctx.params)
File "/usr/local/lib/python3.9/site-packages/click/core.py", line 610, in invoke
return callback(*args, **kwargs)
File "/usr/local/lib/python3.9/site-packages/scanpy_scripts/cmd_utils.py", line 45, in cmd
adata = _read_obj(input_obj, input_format=input_format)
File "/usr/local/lib/python3.9/site-packages/scanpy_scripts/cmd_utils.py", line 87, in _read_obj
adata = sc.read(input_obj, **kwargs)
File "/usr/local/lib/python3.9/site-packages/scanpy/readwrite.py", line 112, in read
return _read(
File "/usr/local/lib/python3.9/site-packages/scanpy/readwrite.py", line 713, in _read
return read_h5ad(filename, backed=backed)
File "/usr/local/lib/python3.9/site-packages/anndata/_io/h5ad.py", line 421, in read_h5ad
d[k] = read_attribute(f[k])
File "/usr/local/lib/python3.9/functools.py", line 877, in wrapper
return dispatch(args[0].__class__)(*args, **kw)
File "/usr/local/lib/python3.9/site-packages/anndata/_io/utils.py", line 183, in func_wrapper
raise AnnDataReadError(
anndata._io.utils.AnnDataReadError: Above error raised while reading key '/layers' of type <class 'h5py._hl.group.Group'> from /.
```

This is likely a Tool Version error. If you use a newer version of a tool with an AnnData object, and then try and use an older version of the tool or other tool in the same toolsuite (Scanpy) later, this will often fail with the above error message. The Scanpy toolsuite is not 'backwards compatable' - few toolsuites are. If this happened while performing a tutorial, we recommend Tutorial Mode as this embeds the correct tool version in each tool button.

{% snippet faqs/galaxy/tutorial_mode.md %}

To fix this in your current history, try re-running the tool with the newer tool version. Or, re-run the prior dataset with an older version.

{% snippet faqs/galaxy/tools_change_version.md %}
