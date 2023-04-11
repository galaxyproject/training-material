---
layout: tutorial_hands_on

title: "Deploying a Beacon v1 in Galaxy"
zenodo_link: ""
questions:
    - What is a Beacon?
    - How do I deploy it?
objectives:
    - Deploy a Beacon
time_estimation: "30m"
key_points:
    - idk bru
contributions:
    authorship:
    - hexylena
    funding:
    - CINECA-Project
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
draft: true
tags:
  - ga4gh
  - beacon
  - git-gat
---

This tutorial will guide you through setting up a Beacon!
TODO: write some more things about beacon

> <warning-title>Beacon v1</warning-title>
> This deploys an older Beacon v1 which was a simpler system.
> The Beacon v1 is more or less deprecated, with users being pushed to Beacon v2 which gives much richer answers, and offers better querying syntax.
{: .warning}

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

## Installing and Configuring

> <hands-on-title>Setting up a Beacon with Ansible</hands-on-title>
>
> 1. Setup the hosts
>
>    {% raw %}
>    ```diff
>    --- a/hosts
>    +++ b/hosts
>    @@ -6,3 +6,12 @@ galaxyservers
>     gat-0.au.training.galaxyproject.eu ansible_user=ubuntu
>     [monitoring]
>     gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +
>    +[beacon]
>    +[beacon:children]
>    +beacon-import
>    +beacon-server
>    +[beacon-server]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +[beacon-import]
>    +gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add hosts"}
>
>    TODO: tip box about children. We get to learn ansible features!
>
> 1. Setup the requirements
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -36,3 +36,7 @@
>       version: 2.1.3
>     - src: galaxyproject.proftpd
>       version: 0.3.1
>    +- name: paprikant.beacon
>    +  src: https://github.com/Paprikant/ansible-role-beacon
>    +- name: paprikant.beacon-importer
>    +  src: https://github.com/Paprikant/ansible-role-beacon_importer
>    {% endraw %}
>    ```
>    {: data-commit="Add requirements"}
>
> 2. Install the role with:
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-galaxy install -p roles -r requirements.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in}
>
> 3. Create the vars file
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/group_vars/beacon.yml
>    @@ -0,0 +1,31 @@
>    +---
>    +postgres_data_dir: /data/beacon/postgresql/data
>    +postgres_init_dir: /data/beacon/postgresql/init
>    +bp_external_binding: 5050 # The default
>    +postgres_user: "{{ beacon_db_user }}"
>    +postgres_pass: "{{ beacon_db_password }}"
>    +postgres_external_binding: "{{ beacon_db_port }}"
>    +# Database Configuration
>    +beacon_db_user: beacon
>    +beacon_db_host: "{{ groups['beacon-server'][0] }}"
>    +beacon_db_password: "{{ beacon_db_password }}"
>    +beacon_db_port: 8080
>    +#galaxy_api_key: This we will set in secrets.
>    +# Information about your beacon (consider filling this out.
>    +beacon_info_title: GA4GH Beacon
>    +beacon_info_beacon_id: your.galaxy.beacon
>    +beacon_info_description: Beacon service hosting datasets from all over the Galaxy
>    +beacon_info_url: https://{{ ansible_inventory_name }}/beacon/
>    +beacon_info_service_group: galaxy-eu
>    +beacon_info_org_id: usegalaxy.aq
>    +beacon_info_org_name: Some Galaxy
>    +beacon_info_org_description: Galaxy community
>    +beacon_info_org_address: 123 Main Street, ZA
>    +beacon_info_org_welcome_url: https://galaxyproject.org/
>    +beacon_info_org_contact_url: https://galaxyproject.org/
>    +beacon_info_org_logo_url: https://galaxyproject.org/images/galaxy-logos/galaxy_project_logo_square.png
>    +beacon_info_org_info: More information about the organisation than just the description can go here.
>    +# Script Configuration
>    +galaxy_api_url: "{{ galaxy_server_hostname }}"
>    +script_user: beacon
>    +script-dir: /home/beacon/script
>    {% endraw %}
>    ```
>    {: data-commit="Add relevant group variables"}
>
>    TODO: tip about 'groups' variable
>
> 3. Add the beacon-server playbook
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/beacon-server.yml
>    @@ -0,0 +1,9 @@
>    +---
>    +- name: Beacon Server
>    +  hosts: beacon-server
>    +  become: true
>    +  become_user: root
>    +  vars_files:
>    +    - group_vars/secret.yml
>    +  roles:
>    +    - paprikant.beacon
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon server playbook"}
>
> 5. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook beacon-server.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
>    TODO: Check that it works
>
{: .hands_on}

## Setting up the Importer

Now that our beacon is running, we need to get data from Galaxy to the Beacon

> <hands-on-title>Setting up the Beacon Importer</hands-on-title>
>
> 1. Add the beacon-import playbook
>
>    {% raw %}
>    ```diff
>    --- /dev/null
>    +++ b/beacon-import.yml
>    @@ -0,0 +1,9 @@
>    +---
>    +- name: Beacon Importer
>    +  hosts: beacon-import
>    +  become: true
>    +  become_user: root
>    +  vars_files:
>    +    - group_vars/secret.yml
>    +  roles:
>    +    - paprikant.beacon-importer
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon importer playbook"}
>
> 1. Edit your `group_vars/secret.yml` and define some random passwords:
>
>    - The API key for your account, which must be an admin
>    
>    {% snippet faqs/galaxy/preferences_admin_api_key.md admin=true %}
>
>    > <code-in-title>Bash</code-in-title>
>    > ```
>    > ansible-vault edit group_vars/secret.yml
>    > ```
>    {: .code-in}
>
>    ```yaml
>    galaxy_api_key: your api key from galaxy
>    ```
>
>    <!-- Ignore this, just for the gat-automation. Vaults are ugly to work with :(
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,13 +1,14 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -62346261323266656232393034396134316636376533376139666437363535393562663838613938
>    -6336666266633563346337623265353935646361326337610a393834333233313461346439376438
>    -63383338346530656561636631666134373238366364363164313166346461383736613162653237
>    -3461363334323431370a656132303965653262386130353332623937376261396530393761353834
>    -38336565666437666436643163363831633331333766653266356163613138393734656465323634
>    -39366362383433366437353534663134313330316337393335383962613961386665633261616237
>    -35366635373063313631323939396164336330356361393464326636353037336461323531336434
>    -35613933303333623031353936393265636130363335376533393335663266313863376135383338
>    -36613464373231623938373434306266373234633036343636633963353361356631363533353066
>    -39323064336237646432323530313065303331326636353334343862373330313133326363363063
>    -38383564636161396435666164643334656435393533643163393434623434656238633631633939
>    -33353232666432376661
>    +34653138616237303566323134333934633262663531653733393332353263663736373737393136
>    +3861366264656366333663303138353266303139356332650a343562633839366135366331666138
>    +38386166346138646364626539383336633763303633623634326235643335636232613061306439
>    +3134626233663563320a363533393066623165376331643663303233396232316161656436626339
>    +64663838343565376465613962663832663862653264633331316363353664343735363163366361
>    +62363633306233316534353737346563353564386235663634306233663332356162646632383934
>    +33613134313231616432303237343933313863623738363330656631373936343966343832636437
>    +39383666303965333762306131393565313366613261376333343630383234336131386165313230
>    +66653830363335303936616364653538613238376235386539643461376432663835303535666462
>    +66656461613530613137393039376234633235353235613064303435663937376437613461333837
>    +30363431326631323736666563633263623966376138656630616464633834313363366239346637
>    +31333562663962376561353034653662623337376665363731393235633065306332623734643264
>    +35353635383966313261616538346661366365636365313631373230383565333037
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon/galaxy passwords to the vault"}
>
>    -->
>
> 1. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook beacon-importer.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
>
> 1. Add the nginx routes
>
>    {% raw %}
>    ```diff
>    --- a/templates/nginx/galaxy.j2
>    +++ b/templates/nginx/galaxy.j2
>    @@ -95,4 +95,14 @@ server {
>          location /reports/ {
>              proxy_pass http://unix:{{ galaxy_config.gravity.reports.bind }}:/;
>          }
>    +
>    +    location /beacon {
>    +        proxy_pass http://{{ groups['beacon-server'][0] }};
>    +        proxy_http_version 1.1;
>    +        proxy_set_header Upgrade $http_upgrade;
>    +        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
>    +        proxy_set_header X-Real-IP $remote_addr;
>    +        proxy_set_header Connection $connection_upgrade;
>    +        proxy_set_header Host $host;
>    +    }
>     }
>    {% endraw %}
>    ```
>    {: data-commit="Add beacon routes to nginx"}
>
> 1. Run the playbook
>
>    > <code-in-title>Bash</code-in-title>
>    > ```bash
>    > ansible-playbook galaxy.yml
>    > ```
>    > {: data-cmd="true"}
>    {: .code-in }
{: .hands_on}

Congratulations, you've set up a Beacon for Galaxy

## Check that it works

????
