---
title: How can I get my container requiring jobs to run in a container?
area: galaxy
box_type: tip
layout: faq
contributors: [hexylena, thanhleviet]
---

Some tools will only run in a container, i.e. they have a container defined in the 'requirements' section of the tool's XML file.
Galaxy will not refuse to run these tools if the container isn't available or if Galaxy isn't configured use containers. Instead it'll run in the host system and likely fail.

## Job Configuration

You can resolve this by configuring your job conf to have destinations that support containers (or even require them.):

The destination must have `docker_enabled` (Or `singularity_enabled`), and you can consider adding `require_container` to make sure the job will fail if the container isn't available. The `docker_volumes` string will allow you to control which volumes are attached to that container;

In TPV configuration (provided by @gtn:thanhleviet) this would look like:

```yaml
  docker:
    inherits: slurm
    scheduling:
      require:
        - docker
    params:
      docker_enabled: true
      require_container: true

  singularity:
    inherits: slurm
    scheduling:
      require:
        - singularity
    params:
      singularity_enabled: true

  podman:
    inherits: slurm
    scheduling:
      require:
        - podman
    params:
      docker_enabled: true
      require_container: true
      docker_volumes: "$galaxy_root:ro,$tool_directory:ro,$job_directory:ro,$working_directory:z,$default_file_path:z"
      docker_sudo: false
      docker_cmd: /usr/bin/podman
      docker_run_extra_arguments: "--userns=keep-id"
```

Or in XML:

```xml
<destination id="docker" runner="local">
    <param id="docker_enabled">true</param>
    <param id="require_container">true</param>
</destionation>

<destination id="podman" runner="local">
    <param id="docker_enabled">true</param>
    <param id="require_container">true</param>
    <param id="docker_sudo">false</param>
    <param id="docker_cmd">/usr/bin/podman</param>
    <param id="docker_run_extra_arguments">--userns=keep-id</param>
    <!-- This will not work until https://github.com/galaxyproject/galaxy/pull/18998 is merged for SELinux users. For now you may want to patch it manually. -->
    <!-- <param id="docker_volumes">$galaxy_root:ro,$tool_directory:ro,$job_directory:ro,$working_directory:z,$default_file_path:z</param> -->
</destination>

<destination id="singularity" runner="local">
    <param id="singularity_enabled">true</param>
    <param id="require_container">true</param>
</destionation>
```


Configuring a tool to use this destination would look like:

```yaml
toolshed.g2.bx.psu.edu/repos/thanhlv/metaphlan4/metaphlan4/4.0.3:
    cores: 12
    mem: cores * 8
    params:
      singularity_enabled: true
```

Or in XML:

```xml
<tools>
    <tool id="toolshed.g2.bx.psu.edu/repos/thanhlv/metaphlan4/metaphlan4/4.0.3" destination="docker"/>
</tools>
```

## Container Resolvers Configuration

If you're using the default `container_resolvers_conf.yml` then there is nothing you need to do. Otherwise you may want to ensure that you have items in there such as `explicit` and `explicit_singularity` among others. See the [galaxy documentation](https://docs.galaxyproject.org/en/master/admin/container_resolvers.html) on the topic.

## Testing

Here is an example of a tool that requires a container, that you can use to test your container configuration:

```xml
<tool name="container-test" id="container" version="5.0" profile="21.09">
        <requirements>
                <container type="docker">ubuntu:22.04</container>
        </requirements>
        <command><![CDATA[
                pwd >> '$output';
                hostname -f >> '$output';
                echo "" >> '$output';
                cat /etc/os-release >> '$output';
                echo "" >> '$output';
                env | sort >> '$output';
]]></command>
        <inputs>
        </inputs>
        <outputs>
                <data name="output" format="txt" label="log" />
        </outputs>
        <help><![CDATA[]]></help>
</tool>
```
