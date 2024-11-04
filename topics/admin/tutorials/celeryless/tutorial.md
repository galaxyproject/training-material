---
layout: tutorial_hands_on

title: "Alternative Celery Deployment for Galaxy"
zenodo_link: ""
questions:
  - What is *required* for Celery to work in Galaxy?
objectives:
  - Setup the bare minimum configuration to get tasks working
  - Avoid deploying, securing, and managing RabbitMQ and Redis and Flower
time_estimation: "1h"
key_points:
  - While a combination of RabbitMQ and Redis is perhaps the most production ready, you can use Postgres as a backend for Celery
  - This significantly simplifies operational complexity, and reduces the attack surface of your Galaxy.
contributions:
  authorship:
  - hexylena
requirements:
  - type: "internal"
    topic_name: admin
    tutorials:
      - ansible
      - ansible-galaxy
subtopic: data
tags:
  - ansible
---

Celery is a new component to the Galaxy world (ca 2023) and is a distributed task queue that *can* be used to run tasks asynchronously. It isn't mandatory, but you might find some features you expect to use to be missing without it.

If you are running a large production deployment you probably want to follow the [Celery+Redis+Flower Tutorial]({% link topics/admin/tutorials/celery/tutorial.md %}).

However, if you are running a smaller Galaxy you may not want to manage deploying Celery (past what Gravity does for you automatically), you may not want to add Redis to your stack, and you may not have need of Flower!

> <agenda-title></agenda-title>
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Configuring Galaxy to use Postgres

AMQP is a message queue protocol which processes can use to pass messages between each other. While a real message queue like RabbitMQ is perhaps the most robust choice, there is an easier option: Postgres

Add the following to your Galaxy configuration to use Postgres:

```bash
amqp_internal_connection: "sqlalchemy+postgresql:///galaxy?host=/var/run/postgresql"
```

# Configuring Celery to use Postgres

Celery would prefer you use Redis (a Key-Value store) as a backend to store results. But we have a database! So let's try using that instead:

```
enable_celery_tasks: true
celery_conf:
  broker_url: null  # This should default to using amqp_internal_connection
  result_backend: "db+postgresql:///galaxy?host=/var/run/postgresql"
  task_routes:
    galaxy.fetch_data: galaxy.external
    galaxy.set_job_metadata: galaxy.external
```

With that we should now be able to [use useful features like](https://docs.galaxyproject.org/en/master/admin/production.html#use-celery-for-asynchronous-tasks):

- Changing the datatype of a collection.
- Exporting histories
- other things!
