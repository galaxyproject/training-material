:whale: Galaxy Docker repository for Exome sequencing training material
====

# Requirements

- [`docker`](https://docs.docker.com/installation/)

# Usage

## Standard usage

To launch:

```
docker run -d -p 8080:80 bgruening/galaxy-exome-seq-training
```

Explanation of this command line:

- `docker run` will run the Image/Container for you.
  In case you do not have the Container stored locally, docker will download it
  for you.
- `-p 8080:80` will make the port 80 (inside of the container) available on
port 8080 on your host.
  Inside the container a Apache Webserver is running on port 80 and that port
  can be bound to a local port on your host computer. With this parameter you
  can access your Galaxy instance via [http://localhost:8080](http://localhost:8080)
- `bgruening/galaxy-exome-seq-training` is the Image/Container name,
that directs docker to the correct path in the [docker index](https://index.docker.io/u/bgruening/galaxy-exome-seq-training/).
- `-d` will start the docker container in daemon mode

For a more detailed description, please consult
the [docker manual](http://docs.docker.io/). It's really worth reading.

Docker images are "read-only", all your changes inside one session will be lost
after restart. This mode is useful to present Galaxy to your colleagues or to
run workshops with it.

To install Tool Shed repositories or to save your data, you need to export the
computed data to the host computer with

```
docker run -d -p 8080:80 -v /home/user/galaxy_storage/:/export/ bgruening/galaxy-exome-seq-training
```

With the additional `-v /home/user/galaxy_storage/:/export/` parameter, Docker
will mount the folder `/home/user/galaxy_storage` into the Container under
`/export/`. A `startup.sh` script, that is usually starting Apache, PostgreSQL
and Galaxy, will recognize the export directory with one of the following outcomes:

- In case of an empty `/export/` directory, it will move the
[PostgreSQL](http://www.postgresql.org/) database, the Galaxy database directory,
Shed Tools and Tool Dependencies and various config scripts to `export/` and
symlink back to the original location.
- In case of a non-empty `/export/`, for example if you continue a previous
session within the same folder, nothing will be moved, but the symlinks will be
created.

This enables you to have different export folders for different sessions - means
real separation of your different projects.

## Stop the Docker container

When you are done with the current Docker container, you can stop it:

- Get the container id or name by executing `docker ps`
- Stop the container with `docker stop <container_id_or_name>`
- (Optional) Delete the container with `docker rm <container_id_or_name>`

## Usage with an interactive session

For an interactive session, execute:

```
docker run -i -t -p 8080:80 bgruening/galaxy-exome-seq-training
```

Then, run the `startup` script by your own, to start PostgreSQL, Apache and Galaxy.

## Launch the Docker container outside the available Docker image

If you want to use directly this repository content outside the available Docker
image:

- Clone this repository
- Move to the current subdirectory (`cd Exome-Seq/docker`)
- Build an image of the docker content with `docker build -t galaxy-exome-seq-training .`
- Run your container with `docker run -p 8080:80 galaxy-exome-seq-training`
- Access your Galaxy instance via [http://localhost:8080](http://localhost:8080)

Once you are done with your container, you can stop and deleted it as described
before.

## Users & Passwords

The Galaxy Admin User has the username `admin@galaxy.org` and the password `admin`.
The PostgreSQL username is `galaxy`, the password is `galaxy` and the database
name is `galaxy` (I know I was really creative ;)).
If you want to create new users, please make sure to use the `/export/` volume.
Otherwise your user will be removed after your docker session is finished.
