---
title: "Enhancing Scientific Training: The Galaxy Training Network's Role in the ELIXIR Training Life-Cycle"
contributions:
  authorship: [bebatut, lisanna]
tags: [gtn]
layout: news
cover: "news/images/2025-07-01-gtn-elixir-life-cycle_training_lifecycle.png"
coveralt: "The image illustrates the training life cycle through a circular diagram featuring five interconnected stages. Starting with the Plan stage, depicted by a checklist icon, it involves identifying and assessing training needs to ensure they align with desired outcomes. The Design stage, represented by a pencil and ruler icon, focuses on creating effective training methods and determining the best delivery approach. In the Develop stage, shown with a computer screen displaying code, training materials such as job aids, activities, presentations, and assessments are created. The Deliver stage, illustrated with a document and folder icon, covers logistical aspects of training delivery, including registration, advertisement, and course execution. Finally, the Evaluate stage, represented by a checklist with checkmarks, involves assessing the training's effectiveness and integrating findings for continuous improvement. The circular arrangement of these stages emphasizes the continuous and iterative nature of the training life cycle."

---

In the rapidly evolving landscape of data science, **continuous learning** and **skill development** are crucial. The Galaxy Training Network (GTN) plays a pivotal role in this educational ecosystem, particularly within the **ELIXIR Training Life-Cycle**. This blog post explores how the GTN contributes to each phase of the life-cycle and aligns with the SPLASH recommendations, ensuring high-quality training for researchers worldwide.

# The Training Life-Cycle

Training in the data sciences often follows a **structured path** from identifying needs to evaluating impact, known as the **Training Life-Cycle**. This cycle helps **training providers** to plan, design, develop, deliver, and evaluate training activities effectively in order to create **high-quality, sustainable training**.

![The image illustrates the training life cycle through a circular diagram featuring five interconnected stages. Starting with the "Plan" stage, depicted by a checklist icon, it involves identifying and assessing training needs to ensure they align with desired outcomes. The "Design" stage, represented by a pencil and ruler icon, focuses on creating effective training methods and determining the best delivery approach. In the "Develop" stage, shown with a computer screen displaying code, training materials such as job aids, activities, presentations, and assessments are created. The "Deliver" stage, illustrated with a document and folder icon, covers logistical aspects of training delivery, including registration, advertisement, and course execution. Finally, the "Evaluate" stage, represented by a checklist with checkmarks, involves assessing the training's effectiveness and integrating findings for continuous improvement. The circular arrangement of these stages emphasizes the continuous and iterative nature of the training life cycle.]({% link news/images/2025-07-01-gtn-elixir-life-cycle_training_lifecycle.png %})

# ELIXIR Training Life-Cycle and SPLASH

To guide trainers through the Training Life-Cycle, the [**ELIXIR Training Platform**](https://elixir-europe.org/platforms/training) developed and maintains the [**SPLASH framework**](https://elixir-europe-training.github.io/ELIXIR-Training-SPLASH/) (Skills: Developing essential skills for trainers and learners; Professional Development: Supporting ongoing professional growth; Learning Assessment: Providing tools and methods for assessing learning outcomes; Support: Offering resources and assistance throughout the training life cycle; Help: Ensuring trainers have access to necessary help and guidance). 

This framework provides support for the different steps of the Training Life-Cycle by offering **recommendations** and linking to **services/tools** supported by ELIXIR.

# Comparison of ELIXIR Services/Tools and GTN Implementations

Let's examine the different steps of the Training Life-Cycle and their associated SPLASH Recommendations, comparing ELIXIR Services/Tools and GTN Implementations.

## Plan

In this stage, **training needs** are identified and assessed to determine whether training is needed and what needs to be trained. This stage helps ensure that the training results in the desired and intended outcomes.

SPLASH Recommendation | ELIXIR Services/Tools | GTN Implementations
--- | --- | ---
**Training Needs Analysis**: Identify specific training gaps and assess requirements. | No dedicated ELIXIR tool exists specifically for needs analysis; often done via surveys at Node level or community consultations. | No dedicated GTN mechanism for collecting training needs. Community contributions driven by perceived gaps but not through a formal system.
**Define Target Audience**: Specify the group benefiting from the training. | Audience (e.g., wet lab researchers, developers) often defined for ELIXIR training events. Some guidance via the ELIXIR Training Platform. | Metadata for audience level (beginner, intermediate, advanced), domain, and prerequisites in **GTN tutorials**.
**Set Learning Objectives**: Outline desired outcomes of the training. | Included in the [**ELIXIR Lesson Template**](https://github.com/elixir-europe-training/ELIXIR-TrP-LessonTemplate-MkDocs). Encouraged by **Train-the-Trainer** events. | Dedicated section for Learning Objectives, aligned with Bloom's taxonomy, in all GTN tutorials.
**Curriculum Development**: Create a detailed curriculum aligned with objectives. | [ELIXIR Learning Paths](https://elixir-europe.org/focus-groups/learning-paths) guiding curriculum structure. | Curated sets of tutorials grouped by topic or goal in [**GTN Learning Pathways**]({% link learning-pathways/index.md %}).
**Train-the-Trainer Programs**: Offer continuous training for trainers. | Well established [**ELIXIR Train-the-Trainer (TtT)**](https://elixir-europe.org/platforms/training/train-the-trainer) to promote pedagogical best practices. | Training and Contributing resources and community onboarding, including materials to support new trainers, including a [**Train the Trainer Learning Pathway**]({% link learning-pathways/train-the-trainers.md %}).

## Design

After identifying the training needs and determining the objectives, an effective **training method** can be designed. This stage includes determining the best **method of delivery**, which may include instructor-led, online, self-study, or a blended method.

SPLASH Recommendation | ELIXIR Services/Tools | GTN Implementations
--- | --- | ---
**Define Learning Outcomes**: Ensure outcomes are SMART and aligned with Bloom's taxonomy. | Promoted in ELIXIR Train-the-Trainer (TtT) and **ELIXIR Lesson Templates**. | Clearly defined Learning Objectives, aligned with Bloom's taxonomy, required in all GTN tutorial.
**Develop Lesson Plans**: Create structured outlines detailing sessions, objectives, and materials. | Encouraged via ELIXIR Course Design Guidelines and Templates. | Consistent template for: GTN tutorials, including structured objectives, content blocks, and teaching aids; [**GTN events**]({% link events/index.md %}); and GTN Learning Pathways
**Select Delivery Format**: Choose the appropriate mode: classroom, online, blended. | Support for various delivery formats across nodes, but no tooling to support format transformations. | Support for **synchronous, asynchronous, and blended learnin**g with **video library** with recorded and subtitled tutorials and **automated tooling** to generate videos from slides, embed interactive components across formats — [Video Toolkit]({% link news/2024/06/14/gtn-video-library.md %}), generate Notebooks, etc
**Design Assessments**: Create quizzes, reflections, or assignments to measure learning outcomes. | [**ELIXIR Training Metrics Database**](https://tmd.elixir-europe.org/world-map) collecting assessment usage data without implementation tools. | Direct integration and **peer-review** of assessments in tutorial structure: **Hands-on exercises**, **Inline question boxes**, **Quizzes and formative checks**.

## Develop

In this stage, materials including **tutorials, activities, presentations, and assessments** are developed. The **training material assets** are created and stored on a training material **repository** or an **e-learning platform**.

SPLASH Recommendation | ELIXIR Services/Tools | GTN Implementations
--- | --- | ---
**Assign Roles and Responsibilities**: Collaborate effectively during material creation. | Collaboration encouraged, particularly through Train-the-Trainer (TtT) and community-driven development. **GitHub-based development** encouraged by the lesson template | GitHub-based development with clear contributor roles (author, editor, reviewer, tester, infrastructure provider, translator, funder) plus automatic testing and content validation. **Editorial board** for topics and learning pathways. [**Individual Contributor / Organizer / Funding page** highlighting contributions]({% link hall-of-fame.md %}).
**Create Diverse Learning Materials**: Use various media types to support different learning styles. | Diverse training materials hosted on TeSS, but it does not offer tools for creating them. | Tutorials combining text and hands-on activities. Support for videos and slide decks. Multimedia integrated using standardized formats and auto-generated content pipelines.
**Ensure FAIR Principles**: Apply Findability, Accessibility, Interoperability, and Reusability. | Extensive guidance and examples in the [FAIR Training Handbook](https://elixir-europe-training.github.io/ELIXIR-TrP-FAIR-training-handbook/#) | [**Perfect FAIRChecker score**]({% link news/2024/05/13/fair.md %}) and follows FAIR-by-design practices via metadata following **Bioschemas**, Zenodo DOI releases, and **ELIXIR TeSS** registration.
**Implement Quality Assurance**: Review, test, and continuously improve training materials. | Peer review encouraged, especially through TtT and community standards, but no formal infrastructure. | Dedicated Editorial Board for topics and learning pathways to oversee content review, style, accuracy, and metadata compliance. Tutorials regularly updated as they are used. Regular updates implemented via GitHub workflows.

## Deliver

This stage focuses on the **delivery of the training**, including logistical aspects like registration, advertisement, and the actual execution of the course, either synchronously or asynchronously if it is an e-learning course.

SPLASH Recommendation | ELIXIR Services/Tools | GTN Implementations
--- | --- | ---
**Prepare for Delivery Format**: Select the appropriate training format and setup (e.g., classroom, online, hybrid). | Recommendations and best practices through TtT and course design materials. | Tutorials designed by default to support:  Self-paced learning, Live virtual or in-person events, Hybrid formats. **Templates for Event Pages**. Tutorials with **delivery recommendations**.
**Ensure Accessibility and Inclusivity**: Accommodate diverse learner needs (e.g., impairments, backgrounds). | Inclusive practices encouraged. Accessibility principles discussed in TtT. Code of Conduct promoted but not centralized. | [Accessibility]({% link accessibility.md %}) (e.g., alt-text, keyboard navigation) embedded. Materials in [multiple formats/languages]({% link faqs/gtn/translations.md %}). [Galaxy Code of Conduct](https://galaxyproject.org/community/coc/) enforced at all events and within the training ecosystem.
**Gather Feedback and Assess Effectiveness**: Collect data during/after the event to improve. | Metrics aggregated in **Training Metrics Database (TMD)** with feedback form templates — but assumes the trainer collects and inputs the data. | [Automated feedback aggregation]({% link feedback.md %}) across tutorials and events, with currently 3,000+ responses, providing real-time insights into learner experience. [Dedicated feedback collection for workshop organizers](https://github.com/galaxyproject/training-material/discussions/1452).
**Support Trainers**: Ensure clear responsibilities and cooperation before/during events. | Guidance available through Train-the-Trainer. Collaboration encouraged, but infrastructure support is decentralized. | Co-teaching and shared authorship encouraged. Pre-configured course handbooks and planning tools for trainers in the Event Pages 
**Adapt to Learners During Delivery**: Respond to questions, pace adaptively, offer support. | Highlighted in TtT workshops and event planning resources. | Structured question prompts, practical guidance, and [**Frequently Asked Questions**]({% link faqs/index.md %}) in each tutorial. **Modular materials** and [**Choose Your Own Tutorial**]({% link news/2022/04/12/cyot.md %}) support, making it easier to adapt on the fly.

## Evaluate

After participants have had a chance to apply the learning, the **effectiveness of the training** should be assessed. These findings are integrated into an overall evaluation and review process to allow improvements to the **training program**.

SPLASH Recommendation | ELIXIR Services/Tools | GTN Implementations
--- | --- | ---
**Impact Assessment**: Evaluate long-term benefits and alignment with broader goals. | [**ELIXIR Impact Toolkit**](https://f1000research.com/documents/12-1127) providing general impact guidance, but not training-specific. Trainers must interpret and adapt it for educational use. | Training-specific impact collected and analyzed via follow-up feedback, tutorial usage statistics, and community engagement metrics.
**Data Collection & Analysis**: Streamline and standardize evaluation data processes. | Templates for standardized collection in the Training Metrics Database (TMD), but assumes trainers manage the process. | **Automated feedback aggregation** from tutorials and events (over 3,000+ responses), [visit (Plausible) analytics](https://plausible.galaxyproject.eu/training.galaxyproject.org?period=12mo&page=/training-material/faqs/), and periodic community reports (e.g. [single-cell]({% link topics/single-cell/community.md %})).
**Report Findings & Improvements**: Communicate outcomes and adjust future trainings. | No centralized or formal reporting mechanism; up to individual trainers or nodes. | Tutorials iteratively updated through GitHub, with visible changelogs, monthly versioning on Zenodo, and community-led issue tracking.

# Conclusion

The Galaxy Training Network is a vital component of the Training Life-Cycle. By adhering to the SPLASH recommendations and utilizing the robust capabilities of the Galaxy platform, the GTN guarantees that researchers have access to high-quality, practical, and sustainable training resources. As the field of data science continues to advance, the GTN's dedication to excellence and collaboration will be essential in empowering the next generation of scientists.

To keep the tables accurate and up-to-date, **please contribute by adding comments to the [following spreadsheets](https://docs.google.com/spreadsheets/d/1WjmjitJfwj3VczKSrVlvPjUeI-FirrpXsjzuC0dKnaU/edit?usp=sharing)!**
