---
layout: tutorial_hands_on

title: "Deploying a Beacon v1 in Galaxy"
zenodo_link: ""
questions:
    - What is a Beacon?
    - How do I deploy it?
    - Is v1 the same as v2?
objectives:
    - Deploy a Beacon
time_estimation: "30m"
key_points:
    - While deprecated, Beacon v1 is easy to deploy
    - It can also tick some boxes for grants!
contributions:
    authorship:
    - hexylena
    editing:
    - shiltemann
    funding:
    - CINECA-Project
subtopic: features
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
tags:
  - ga4gh
  - beacon
  - git-gat
---

This tutorial will guide you through setting up a [GA4GH Beacon](https://beacon-project.io/)!

> The Beacon Project is developed under a Global Alliance for Genomics and Health (GA4GH) Iniciative for the federated discovery of genomic data in biomedical research and clinical applications. 
{: .quote}

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
>    @@ -9,3 +9,12 @@ gat-0.eu.galaxy.training ansible_connection=local ansible_user=ubuntu
>     
>     [sentryservers]
>     gat-0.eu.training.galaxyproject.eu ansible_connection=local ansible_user=ubuntu
>    +
>    +[beacon]
>    +[beacon:children]
>    +beacon_import
>    +beacon_server
>    +[beacon_server]
>    +gat-0.eu.galaxy.training ansible_connection=local ansible_user=ubuntu
>    +[beacon_import]
>    +gat-0.eu.galaxy.training ansible_connection=local ansible_user=ubuntu
>    {% endraw %}
>    ```
>    {: data-commit="Add hosts"}
>
>    > <tip-title>What are these :children?</tip-title>
>    > Here we use some of the more advanced features of Ansible's Inventory system.
>    > We declare a host group called 'beacon' with no hosts of its own.
>    > 
>    > Then we declare that this beacon group has two children: beacon_import, and beacon_server. We can then define host groups for those two entries with as many different hosts as we need. This makes it very easy to scale up configuration.
>    > 
>    > Here we will be using that feature to declare some 'beacon variables', which will be shared between the beacon_server and beacon_importer. Because they're children of 'beacon', they'll inherit any group variables defined for `group_vars/beacon.yml`.
>    {: .tip}
>
> 1. Setup the requirements
>
>    {% raw %}
>    ```diff
>    --- a/requirements.yml
>    +++ b/requirements.yml
>    @@ -60,3 +60,8 @@
>     # Our FTP Server
>     - src: galaxyproject.proftpd
>       version: 0.3.1
>    +# Beacon support
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
>    +beacon_db_host: "{{ groups['beacon_server'][0] }}"
>    +beacon_db_password: "{{ vault_beacon_db_password }}"
>    +beacon_db_port: 9001
>    +#galaxy_api_key: This we will set in secrets.
>    +# Information about your beacon (consider filling this out.
>    +beacon_info_title: GA4GH Beacon
>    +beacon_info_beacon_id: your.galaxy.beacon
>    +beacon_info_description: Beacon service hosting datasets from all over the Galaxy
>    +beacon_info_url: https://{{ groups['beacon_server'][0] }}/beacon/
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
>    +galaxy_api_url: "https://{{ groups['galaxyservers'][0] }}"
>    +script_user: beacon
>    +script_dir: /home/beacon/script
>    {% endraw %}
>    ```
>    {: data-commit="Add relevant group variables"}
>
>    > <tip-title>groups?</tip-title>
>    > Here we again use some advanced features of Ansible's inventory system. Ansible knows the name of every hostname in the inventory. Now that we want to point the beacon configuration, either at the database which should be on `beacon_server`, or at Galaxy in `galaxyservers`, we ask the `groups` variable for what the inventory looks like. We use `[0]` to pull out the first hostname we find for both of those groups.
>    {: .tip}
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
>    +  hosts: beacon_server
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
>    +  hosts: beacon_import
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
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,13 +1,15 @@
>    +galaxy_api_key: your-galaxy-api-key
>    +vault_beacon_db_password: some-super-secret-password
>    ```
>
>    <!-- Ignore this, just for the gat-automation. Vaults are ugly to work with :(
>
>    {% raw %}
>    ```diff
>    --- a/group_vars/secret.yml
>    +++ b/group_vars/secret.yml
>    @@ -1,23 +1,26 @@
>     $ANSIBLE_VAULT;1.1;AES256
>    -31383961303162656266343866376462343961326163356238636565386630363164353330643862
>    -3631663631383136356436363635373536323536353534650a393833343465396339373731663238
>    -38383763326230346561633061356634643566663863383838363765326162386133636330656136
>    -3935643333353365660a393932313965373166313838373739306234306634313139656334343839
>    -38393636396164373266646534326436316333313638366234323564323936643863633533366537
>    -66393535333263383837386538636530303137623761333161383765643938396232633863353838
>    -66363464613466353833663833613339666634333734343163326162393031363530643363396238
>    -36306562333763333435646432613565323537373032643335336639353762646662343865663266
>    -33646334623838333230343861353632626461363436343963636664393965613736306162396235
>    -64373239373739613164623165653832653139646433373531323562373534393563363962373634
>    -30393161303439376138653066646230303464353336656431663137366630323830316561356137
>    -65363663663165323666343335626238666334346664373938333839323161353131643336343738
>    -31646330386239363165643866316431386631643139613065663835323036353631646365333530
>    -37306161656562316665643666363132633766386330383864663937316233376261343264633133
>    -39356438623533376539393662326364666434386330626437646235386261323632396661656661
>    -31613232316237626562646330626564333566633266653762356338646236353137343863633634
>    -31623563393933656635323165326638376634383730393131373963353262366337383632666666
>    -64643264333463373334343433373564326632393637626237353336373837373438303231636633
>    -31666636316533313336613532643437646362636431306262623361326532613238303732396266
>    -62346131663135313935393439343433333439623335623131653066316561646631616361303632
>    -66633935356333643539323461333734376466316262656330396662346232333131313864343636
>    -6634613432653361653133383862383338353432666136303465
>    +35326661346439306337366130643630636230643766333533343364613135326265643638663832
>    +3235653731653438383234386233393939376235666163610a613436623232623036316538356338
>    +39313965393163376666316138353631633537636438336537656466366662393539346634653134
>    +6539316435666632300a383730303262333263363161316533613132653464316536653439323331
>    +61376430336336663031393064396462653130336534623664643362663432613965343237326335
>    +37666237623961663763336362636166666233343438336633313765363539646433383736373433
>    +39303735616131663239313337343637663930373638356564356531323066363131396536333030
>    +64326264633936346533643661336133336633363031316236633065313635313762363461366363
>    +62353061653130366562666633616466646436636361316265313131383234613430656461656466
>    +62393264393532323836333536363062653331393738383861656462656539316661336132663739
>    +31316632376436656639323530616432356131663763333961343764616137376133643061383761
>    +66646461613866653662306239386331396364653132656536623938653134653931363366623231
>    +39626362633136353636383539363966363937616463306466643033333038343162366239336332
>    +30393062613161623535623130616635333330356462623363303231653635613263306635663466
>    +39623433383031393935323362353830346466353961636461386566333230316431623532623361
>    +38636233346335623234396233623662363165616339633536343436626633663236653365653065
>    +37363633393463386233643437633662663135613133616561346361633663393163343636663163
>    +36353333393132656265656432616239356462343439633939303339333633383935346233326664
>    +30303963303365313436346334666630373634323639313635333332656532393965346334363338
>    +63633962383732643063316366636334613236343638353233313962326435616661623230643531
>    +31366230643263393063653634383465303636643333363732336438623063333664326362363561
>    +63343538356536626232623838373138336130663435373830376635366566666532303830633834
>    +66353834356566636333393137643836336130653138633561613530363037393438653363623131
>    +34643433643965643366373335316133323830386635653036613061316231326561346631623865
>    +66343936666137376339313638313864353335343538643763366661333933356566
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
>    @@ -115,4 +115,14 @@ server {
>     
>     	{{ tiaas_nginx_routes }}
>     
>    +	location /beacon {
>    +		proxy_pass http://{{ groups['beacon_server'][0] }}:5050;
>    +		proxy_http_version 1.1;
>    +		proxy_set_header Upgrade $http_upgrade;
>    +		proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
>    +		proxy_set_header X-Real-IP $remote_addr;
>    +		proxy_set_header Connection "upgrade";
>    +		proxy_set_header Host $host;
>    +	}
>    +
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

Congratulations, you've set up a Beacon v1 for Galaxy! Go check it out at `/beacon/` on your server.

{% snippet topics/admin/faqs/missed-something.md step=17 %}

{% snippet topics/admin/faqs/git-gat-path.md tutorial="beacon" %}
