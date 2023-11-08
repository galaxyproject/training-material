---
layout: tutorial_hands_on

title: "Galaxy Monitoring with Reports"
zenodo_link: ""
questions:
  - How to monitor a Galaxy service with the Reports application?
objectives:
  - Setup and start the Galaxy reports app.
time_estimation: "30m"
key_points:
  - Galaxy supports pluggable monitoring extensions.
  - The Reports webapp is one option to monitor your system.
contributors:
  - natefoo
  - bgruening
  - slugger70
  - hexylena
subtopic: monitoring
tags:
  - ansible
  - monitoring
  - git-gat
  - broken
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
---

> <warning-title>Currently Broken, Requires Separate Domain</warning-title>
> Reports does not work, under a path prefix (the default setup that most
> people will use.) It is completely broken and the developers have no plans to fix it in the near term.
> See
> [galaxyproject/galaxy#15966](https://github.com/galaxyproject/galaxy/issues/15966) for more details.
>
> However, it should still function with a separate domain, if that is possible
> for your setup. Otherwise, it **will not work.** If you wish to follow this
> tutorial, please be aware of this.
{: .warning}

The reports application gives some pre-configured analytics screens. These are very easy to setup and can help with debugging issues in Galaxy.

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="reports" %}

# Setting up Reports

The reports application is included with the Galaxy codebase and this tutorial assumes you've already done all of the setup required for Galaxy, systemd, uWSGI, and NGINX.

> <hands-on-title>Setup Reports</hands-on-title>
>
>
> 1. First we add a basic configuration of the Reports app to the playbook templates. Create `templates/galaxy/config/` folder, if it doesn't exist, and create `templates/galaxy/config/reports.yml` with the following contents:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/galaxy/config/reports.yml
>    @@ -0,0 +1,4 @@
>    +reports:
>    +    database_connection: "{{ galaxy_config.galaxy.database_connection }}"
>    +    file_path: "{{ galaxy_config.galaxy.file_path }}"
>    +    template_cache_path: "{{ galaxy_mutable_data_dir }}/compiled_templates/reports/"
>    {% endraw %}
>    ```
>    {: data-commit="Setup reports config file"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. In your `galaxyservers` group variables file, tell the playbook to deploy the reports configuration file, and gravity to manage reports:
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -148,6 +148,11 @@ galaxy_config:
>             pools:
>               - job-handlers
>               - workflow-schedulers
>    +    reports:
>    +      enable: true
>    +      url_prefix: /reports
>    +      bind: "unix:{{ galaxy_mutable_config_dir }}/reports.sock"
>    +      config_file: "{{ galaxy_config_dir }}/reports.yml"
>     
>     galaxy_job_config_file: "{{ galaxy_config_dir }}/galaxy.yml"
>     
>    @@ -168,6 +173,8 @@ galaxy_config_templates:
>         dest: "{{ galaxy_config.galaxy.dependency_resolvers_config_file }}"
>       - src: templates/galaxy/config/job_resource_params_conf.xml.j2
>         dest: "{{ galaxy_config.galaxy.job_resource_params_file }}"
>    +  - src: templates/galaxy/config/reports.yml
>    +    dest: "{{ galaxy_config.gravity.reports.config_file }}"
>     
>     galaxy_extra_dirs:
>       - /data
>    {% endraw %}
>    ```
>    {: data-commit="Enable gravity to manage reports"}
>
>
> 4. Then we need to tell NGINX it should serve our Reports app under `<server_url>/reports` url. Edit your `templates/nginx/galaxy.j2` file, and within the server block, add a block for proxying the reports application. It should look like:
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -103,4 +103,9 @@ server {
>     		proxy_set_header Upgrade $http_upgrade;
>     		proxy_set_header Connection "upgrade";
>     	}
>    +
>    +	location /reports/ {
>    +		proxy_pass http://{{ galaxy_config.gravity.reports.bind }}:/;
>    +	}
>    +
>     }
>    {% endraw %}
>    ```
>    {: data-commit="Configure a location block for Reports in NGINX"}
>
> 5. Run the playbook:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 6. The reports application should be available, under [`/reports`](https://my.gat.galaxy.training/?path=/reports)
{: .hands_on}

> ```bash
> 1.sh
> ```
> {: data-test="true"}
{: .hidden}

> <comment-title>Insecure!</comment-title>
> But notice that your Reports server is not secured! Check out the [External Authentication]({% link topics/admin/tutorials/external-auth/tutorial.md %}) tutorial for information on securing Reports.
{: .comment}

{% snippet topics/admin/faqs/git-commit.md page=page %}

{% snippet topics/admin/faqs/missed-something.md step=15 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="reports" %}
