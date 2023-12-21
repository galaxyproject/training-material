---
layout: tutorial_hands_on

title: "Monitoring Galaxy and Pulsar with Sentry"
zenodo_link: ""
questions:
objectives:
  - Have an understanding of Sentry
  - Install Sentry
  - Configure Galaxy and Pulsar to send errors to Sentry
  - Monitor performance with Sentry
time_estimation: "1h"
key_points:
contributions:
  authorship:
  - mvdbeek
  editing:
  - hexylena
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
      - pulsar
subtopic: monitoring
tags:
  - ansible
  - git-gat
---

# Overview

Sentry is an error tracking software that helps admins and developers monitor and diagnose issues in their applications. It provides real-time alerts for errors and allows users to capture context information about each error, such as stack traces and user feedback. It is often possible to find and fix errors before users report them. Galaxy and Pulsar can log issues and failing tool runs to Sentry.


> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="sentry" %}

We're going to set up a local Sentry instance using docker-compose and connect Galaxy and Pulsar to that Sentry instance. Alternatively, you can use the hosted Sentry at https://sentry.io/.

# Installing and Configuring

To proceed from here it is expected that:

> <comment-title>Requirements for Running This Tutorial</comment-title>
>
> 1. You have set up a working Galaxy instance as described in the [ansible-galaxy](../ansible-galaxy/tutorial.html) tutorial.
>
{: .comment}

# Installing and Configuring

First we need to add our new Ansible role to `requirements.yml`:

> <hands-on-title>Set up Sentry with Ansible</hands-on-title>
>
> 1. In your working directory, add the roles to your `requirements.yml`
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -54,3 +54,6 @@
>     # Training Infrastructure as a Service
>     - src: galaxyproject.tiaas2
>       version: 2.1.5
>    +# Sentry
>    +- name: mvdbeek.sentry_selfhosted
>    +  src: https://github.com/mvdbeek/ansible-role-sentry/archive/main.tar.gz
>    {% endraw %}
>    ```
>    {: data-commit="Add requirement" data-ref="add-req"}
>
>    {% snippet topics/admin/faqs/diffs.md %}
>
> 2. Install the roles with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true" data-ref="req-install"}
>    {: .code-in}
>
> 3. Create a new playbook, `sentry.yml` with the following:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/sentry.yml
>    @@ -0,0 +1,7 @@
>    +- hosts: sentryservers
>    +  become: true
>    +  pre_tasks:
>    +    - pip:
>    +        name: docker-compose
>    +  roles:
>    +    - mvdbeek.sentry_selfhosted
>    {% endraw %}
>    ```
>    {: data-commit="Setup the sentry playbook"}
>
>    During this tutorial we will install everything on the same host, but often one keeps the monitoring infrastructure (Grafana, InfluxDB, Sentry) on a separate host.
>
> 4. Edit the inventory file (`hosts`) an add a group for Sentry like:
>
>    {% raw %}
>    ```diff
>    --- a/hosts
>    +++ b/hosts
>    @@ -6,3 +6,6 @@ galaxyservers
>     gat-0.oz.galaxy.training ansible_user=ubuntu
>     [monitoring]
>     gat-0.eu.galaxy.training ansible_connection=local ansible_user=ubuntu
>    +
>    +[sentryservers]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add the monitoring host"}
>
>    **Ensure that the hostname is the full hostname of your machine.**
>
>    Sentry requires its own (sub)domain. For the admin training we have set up the sentry.gat-N.eu.galaxy.training subdomain. If you run this tutorial outside of the training and you cannot obtain a domain or subdomain for sentry you can use the free [Duck DNS](https://www.duckdns.org/) service to map an IP address to a domain name.
>
> 5. Edit the file `group_vars/sentryservers.yml` and set the following variables:
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/sentryservers.yml
>    @@ -0,0 +1,6 @@
>    +sentry_version: 23.3.1
>    +sentry_url: "https://{{ sentry_domain }}"
>    +sentry_docker_compose_project_folder: /srv/sentry
>    +sentry_superusers:
>    +  - email:  admin@example.com
>    +    password: "{{ vault_sentry_password }}"
>    {% endraw %}
>    ```
>    {: data-commit="Configure Sentry"}
>
> 9. We will add an associated admin password to the vault, do that now:
>
>    ><code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```yaml
>    vault_sentry_password: 'some-super-secret-password'
>    ```
>
>    <!-- Ignore this, just for the gat-automation. Vaults are ugly to work with :(
>
>    +vault_sentry_password: password
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,22 +1,23 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -66653366633665383231303635643739396466653465626163633662623230643534346666303434
>    -3531356339646666623364383435306363373534373833370a333964306130653236373438373264
>    -34366433396135303064643932313135643064373034323865363939333565623238386161333465
>    -3932623162353034370a653663333363383938393936623663343639376430653366646662353830
>    -34383964326161383966306361633162653366353162633564623738353137333936313232363764
>    -38383962396637636136316366633936316365643038333030333932336365326335373737663834
>    -39363639303764343430336435313663313762353335366562656232646334356438343535303032
>    -35636131323163323036383163363964643237313137333131303737346662373233663562343031
>    -39396531306539623661376166313534303462623532323334393736316634616262313633626637
>    -37393066663239303564356232366466333334316565363631363230626665333039386133383933
>    -31656234363064646334623938383637393534313638643635396662643163346633346237353664
>    -34666161373764663166653737323736373261363565666633653831363833393264333666356633
>    -30373062333638366666316132623736386336383933326539663065653538386363323766313664
>    -63323737396264306439303366343834326164393763376663316366353766663131303966393039
>    -61343732313865613761623733356465336263633963376161663931653864383162393839643834
>    -34373239386461666462643162613536316166656135323332316163336635393035653164363138
>    -34336265383033663436306436663536383238366665653239313661666431666339323364646633
>    -64366361653034393962316562373835623631356339313062393933633634653735656636373031
>    -65376162633866623733623961303664663931636130373936326461623465363932373165356563
>    -36323838323865326464346337393164343765613333663830636234666663303535616532323936
>    -3230
>    +31383961303162656266343866376462343961326163356238636565386630363164353330643862
>    +3631663631383136356436363635373536323536353534650a393833343465396339373731663238
>    +38383763326230346561633061356634643566663863383838363765326162386133636330656136
>    +3935643333353365660a393932313965373166313838373739306234306634313139656334343839
>    +38393636396164373266646534326436316333313638366234323564323936643863633533366537
>    +66393535333263383837386538636530303137623761333161383765643938396232633863353838
>    +66363464613466353833663833613339666634333734343163326162393031363530643363396238
>    +36306562333763333435646432613565323537373032643335336639353762646662343865663266
>    +33646334623838333230343861353632626461363436343963636664393965613736306162396235
>    +64373239373739613164623165653832653139646433373531323562373534393563363962373634
>    +30393161303439376138653066646230303464353336656431663137366630323830316561356137
>    +65363663663165323666343335626238666334346664373938333839323161353131643336343738
>    +31646330386239363165643866316431386631643139613065663835323036353631646365333530
>    +37306161656562316665643666363132633766386330383864663937316233376261343264633133
>    +39356438623533376539393662326364666434386330626437646235386261323632396661656661
>    +31613232316237626562646330626564333566633266653762356338646236353137343863633634
>    +31623563393933656635323165326638376634383730393131373963353262366337383632666666
>    +64643264333463373334343433373564326632393637626237353336373837373438303231636633
>    +31666636316533313336613532643437646362636431306262623361326532613238303732396266
>    +62346131663135313935393439343433333439623335623131653066316561646631616361303632
>    +66633935356333643539323461333734376466316262656330396662346232333131313864343636
>    +6634613432653361653133383862383338353432666136303465
>    {% endraw %}
>    ```
>    {: data-commit="Add sentry passwords to the vault"}
>
>    -->
>
> 6. Add the nginx routes
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/templates/nginx/sentry.j2
>    @@ -0,0 +1,20 @@
>    +server {
>    +	# Listen on port 443
>    +	listen        *:443 ssl;
>    +	# The virtualhost is our domain name
>    +	server_name   "{{ sentry_domain }}";
>    +
>    +	# Our log files will go here.
>    +	access_log  syslog:server=unix:/dev/log;
>    +	error_log   syslog:server=unix:/dev/log;
>    +
>    +	location / {
>    +		# This is the backend to send the requests to.
>    +		proxy_pass "http://localhost:9000";
>    +
>    +		proxy_set_header Host $http_host;
>    +		proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
>    +		proxy_set_header X-Forwarded-Proto $scheme;
>    +		proxy_set_header Upgrade $http_upgrade;
>    +	}
>    +}
>    {% endraw %}
>    ```
>    {: data-commit="Add nginx server"}
>
> 6. And make sure the sentry nginx configuration is deployed
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -216,6 +216,7 @@ nginx_servers:
>       - redirect-ssl
>     nginx_ssl_servers:
>       - galaxy
>    +  - sentry
>     nginx_enable_default_server: false
>     nginx_conf_http:
>       client_max_body_size: 1g
>    {% endraw %}
>    ```
>    {: data-commit="Deploy nginx server"}
>
> 7. Run the sentry playbook to deploy sentry and the galaxy playbook to update the nginx configuration.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook sentry.yml galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 8. Generate a project for Galaxy in Sentry
>  Go to the domain you configured for your Sentry instance. You need to log in with the username and admin you've set up in `group_vars/sentryservers.yml`. Click "continue" on the next page. Click "Projects", "Create Project", "Python", select "I'll create my own alerts later", and set "galaxy" as the Project Name. You'll see your project dsn that will look like `https://b0022427ee5345a8ad4cb072c73e62f4@sentry.gat-N.eu.galaxy.training/2`. We will need this string to let Galaxy know where to send data to. To avoid requesting an additional certificate for communication between Galaxy and Sentry we've set up communication via localhost:9000, so you can manually change the @ portion to localhost:9000.
>
> 9. We will add the galaxy project dsn to the vault. Edit your `group_vars/secret.yml` and add the sentry dsn.
>
>    ><code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```yaml
>    vault_galaxy_sentry_dsn: 'https://b0022427ee5345a8ad4cb072c73e62f4@localhost:9000/2'
>    ```
>
> 9. Edit `group_vars/galaxyservers.yml` to reference the new vault secret:
>
>    This will let Galaxy know that captured logs should be sent to our Sentry instance.
>    We will also enable sending performance metrics to Sentry by setting the `sentry_traces_sample_rate` to `0.5`. This will send half of all transactions to Sentry. In a production environment you would reduce this to a smaller percentage of transactions.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -116,6 +116,8 @@ galaxy_config:
>         # Monitoring
>         statsd_host: localhost
>         statsd_influxdb: true
>    +    sentry_dsn: "{{ vault_galaxy_sentry_dsn }}"
>    +    sentry_traces_sample_rate: 0.5
>       gravity:
>         process_manager: systemd
>         galaxy_root: "{{ galaxy_root }}/server"
>    {% endraw %}
>    ```
>    {: data-commit="Configure Galaxy to report to Sentry"}
>
> 10. Run the galaxy playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on }

## Generate an error

Galaxy has a built in route that intentionally generates and error.
Just visit: [/error](https://my.gat.galaxy.training/?path=/error)

> <hands-on-title>Open the Galaxy Project in Sentry</hands-on-title>
> 1. Go to your Sentry instance and click on issues. You should see a couple of issues,
> one them should be the "Fake error" exception we generated by visiting https://galaxy.example.org/error.
{: .hands_on }

## Sending tool error reports to Sentry

In addition to sending logging errors to Sentry you can also collect failing tool runs in Sentry. For this we will set up the error reporting configuration file and reference it in galaxy.yml. The `user_submission` parameter controls whether all reports will be collected in Sentry (when set to `false`) or only those that have been reported manually (when set to `true`). For testing purposes we'll also add a tool that will fail running so we can test that submitting tool errors to Sentry works as expected.

> <hands-on-title>Update Galaxy config to send tool error reports</hands-on-title>
>
> 1. Create the  `files/galaxy/config/error_reports.yml` file.
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/config/error_reports.yml
>    @@ -0,0 +1,2 @@
>    +- type: sentry
>    +  user_submission: false
>    {% endraw %}
>    ```
>    {: data-commit="Configure error reporting"}
>
> 2. Create a testing tool in `files/galaxy/tools/job_properties.xml`.
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/files/galaxy/tools/job_properties.xml
>    @@ -0,0 +1,65 @@
>    +<tool id="job_properties" name="Test Job Properties" version="1.0.0">
>    +    <stdio>
>    +        <exit_code range="127" level="fatal" description="Failing exit code." />
>    +    </stdio>
>    +    <version_command>echo 'v1.1'</version_command>
>    +    <command><![CDATA[
>    +#if $thebool
>    +    echo 'The bool is true' &&
>    +    echo 'The bool is really true' 1>&2 &&
>    +    echo 'This is a line of text.' > '$out_file1' &&
>    +    cp '$out_file1' '$one' &&
>    +    cp '$out_file1' '$two' &&
>    +    sleep $sleepsecs
>    +#else
>    +    echo 'The bool is not true' &&
>    +    echo 'The bool is very not true' 1>&2 &&
>    +    echo 'This is a different line of text.' > '$out_file1' &&
>    +    sleep $sleepsecs &&
>    +    sh -c 'exit 2'
>    +#end if
>    +#if $failbool
>    +    ## use ';' to concatenate commands so that the next one is run independently
>    +    ## of the exit code of the previous one
>    +    ; exit 127
>    +#end if
>    +    ]]></command>
>    +    <inputs>
>    +        <param name="sleepsecs" type="integer" value="0" label="Sleep this many seconds"/>
>    +        <param name="thebool" type="boolean" label="The boolean property" />
>    +        <param name="failbool" type="boolean" label="The failure property" checked="false" />
>    +    </inputs>
>    +    <outputs>
>    +        <data name="out_file1" format="txt" />
>    +        <collection name="list_output" type="list" label="A list output">
>    +            <data name="one" format="txt" />
>    +                <has_line line="The bool is true" />
>    +            </assert_stdout>
>    +            <assert_stderr>
>    +                <has_line line="The bool is really true" />
>    +            </assert_stderr>
>    +            <assert_command_version>
>    +                <has_text text="v1.1" />
>    +            </assert_command_version>
>    +        </test>
>    +        <test expect_exit_code="2">
>    +            <param name="thebool" value="false" />
>    +            <output name="out_file1" file="simple_line_alternative.txt" />
>    +            <assert_command>
>    +                <has_text text="very not" />
>    +            </assert_command>
>    +            <assert_stdout>
>    +                <has_line line="The bool is not true" />
>    +            </assert_stdout>
>    +            <assert_stderr>
>    +                <has_line line="The bool is very not true" />
>    +            </assert_stderr>
>    +        </test>
>    +        <test expect_exit_code="127" expect_failure="true">
>    +            <param name="thebool" value="true" />
>    +            <param name="failbool" value="true" />
>    +        </test>
>    +    </tests>
>    +    <help>
>    +    </help>
>    +</tool>
>    {% endraw %}
>    ```
>    {: data-commit="Configure a tool"}
>
> 3. Edit `group_vars/galaxyservers.yml` to reference the `error_reports.yml` file and the new testing tool.
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/galaxyservers.yml
>    +++ b/group_vars/galaxyservers.yml
>    @@ -118,6 +118,7 @@ galaxy_config:
>         statsd_influxdb: true
>         sentry_dsn: "{{ vault_galaxy_sentry_dsn }}"
>         sentry_traces_sample_rate: 0.5
>    +    error_report_file: "{{ galaxy_config_dir }}/error_reports_file.yml"
>       gravity:
>         process_manager: systemd
>         galaxy_root: "{{ galaxy_root }}/server"
>    @@ -170,6 +171,8 @@ galaxy_config_files:
>         dest: "{{ galaxy_config.galaxy.themes_config_file }}"
>       - src: files/galaxy/config/tpv_rules_local.yml
>         dest: "{{ tpv_mutable_dir }}/tpv_rules_local.yml"
>    +  - src: files/galaxy/config/error_reports.yml
>    +    dest: "{{ galaxy_config.galaxy.error_report_file }}"
>     
>     galaxy_config_templates:
>       - src: templates/galaxy/config/container_resolvers_conf.yml.j2
>    @@ -191,6 +194,7 @@ tpv_privsep: true
>     
>     galaxy_local_tools:
>     - testing.xml
>    +- job_properties.xml
>     
>     # Certbot
>     certbot_auto_renew_hour: "{{ 23 |random(seed=inventory_hostname)  }}"
>    {% endraw %}
>    ```
>    {: data-commit="Deploy files, error reporting"}
>
> 4. Run the galaxy playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on }

## Generating a tool error

To generate a tool error, run the job properties testing tool and set the `failbool` parameter to true.

> <hands-on-title>Open the Galaxy Project in Sentry</hands-on-title>
> 1. Go to your Sentry instance and click on issues. You should see an issue for the tool run error.
{: .hands_on }

## Reporting errors from the Pulsar server

It is also possible to report errors from the Pulsar server. You can either use the Galaxy project we created before in Sentry, or we can create a new project for Pulsar. We recommend creating a separate Pulsar project. Since the Pulsar server runs on a remote VM for this to work you need a valid certificate for the Sentry domain and you cannot use localhost.

> <hands-on-title>Add Sentry connection to Pulsar</hands-on-title>
> 1. Create a new dsn by creating a new pulsar project in Sentry.
> 2. We will add the project dsn to the vault. Edit your `group_vars/secret.yml` and add the sentry dsn.
>
>    ><code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```yaml
>    vault_pulsar_sentry_dsn: 'https://f2a8a00d30224c2c9800a8f79194a32a@{{ groups['sentryservers'][0] }}/3'
>    ```
> 3. Add the sentry dsn to the pulsar group variables.
>    {% raw %}
>    ```diff
>    --- a/group_vars/pulsarservers.yml
>    +++ b/group_vars/pulsarservers.yml
>    @@ -45,6 +45,7 @@ pulsar_yaml_config:
>           - type: conda
>             auto_init: true
>             auto_install: true
>    +  sentry_dsn: "{{ vault_pulsar_sentry_dsn }}"
>     
>     # Pulsar should use the same job metrics plugins as Galaxy. This will automatically set `job_metrics_config_file` in
>     # `pulsar_yaml_config` and create `{{ pulsar_config_dir }}/job_metrics_conf.yml`.
>    {% endraw %}
>    ```
>    {: data-commit="Configure pulsar for error reporting"}
>
> 4. Run the pulsar playbook.
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
{: .hands_on }

Pulsar should now be set up to report errors to Sentry.

{% snippet topics/admin/faqs/missed-something.md step=12 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="sentry" %}
