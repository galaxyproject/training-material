---
layout: doc
title: GTN GenAI Policy
subtitle: |
    The GTN is committed to providing high-quality human-written and human-curated tutorials.
    As such, the use of generative AI tools to create GTN training content is not allowed.
contributions:
    authorship: [shiltemann, hexylena]
    reviewing: [bebatut, bgruening]
---


## Rules

- It is **not permitted** to use generative AI to:
  - create GTN content (tutorials, news posts, events, learning pathways, images, and any other user-facing content)
  - communicate about your contribution (GitHub PR descriptions, comments, issues, reviews, chat messages, etc)

- LLMs **may be used** to answer questions, analyze, distill, refine, check, suggest, review.
  - However, we encourage contributors to avoid them in favour of [existing documentation](https://training.galaxyproject.org/training-material/topics/contributing/) and our [chat forum](https://matrix.to/#/%23Galaxy-Training-Network_Lobby%3Agitter.im).
  - All GTN content and all communication must be in **your own words**.

- For contributions to the **codebase** (e.g. `bin`, `_plugins`) you may be *assisted* by generative AI, under the following conditions:
  - You must **understand and take responsibility** for every line of code in your contribution.
  - You must be able to **answer technical questions in review** and make any requested changes.
  - You must **disclose and describe** the use of generative AI in your PR description.
  - You must write code comments and documentation, yourself. Describing not only what the code does, but why you chose to do it this way.
  - You must write the PR/issue descriptions and any other (GitHub) communications yourself.
  - See the [FAQ](#why-a-technical-codebase-exception) for rationale

- Some parts of the GTN codebase are a mix of technical framework and user-facing content (e.g. _layouts folder).
  - Please use your best judgment to adhere to the spirit of this policy when editing these files; i.e. any time you are teaching or communicating with learners it should be in your own words, technical updates may be AI assisted under the conditions outlined above.
  - If in doubt, please contact us on [Matrix](https://matrix.to/#/%23Galaxy-Training-Network_Lobby%3Agitter.im) to discuss your planned contribution with us

- Every **contributor is responsible** for ensuring that their contribution does not violate our policy or license
- Since AI generated information is frequently misleading or false, we cannot offer support on anything referencing AI output.
- Maintainers are encouraged to deprioritise or reject contributions that they suspect have been made in violation of this policy.


## Motivation

### The GTN is for Humans

Learners already have the ability to ask LLMs to create tutorials for them; The GTN is a resource for human written, curated and verified tutorials.

### We value diverse voices

LLMs [distort our written language](https://arxiv.org/abs/2603.18161), reduce [linguistic diversity](https://www.nature.com/articles/s41562-026-02550-0) and result in [semantic ablation](https://www.theregister.com/software/2026/02/16/semantic-ablation-why-ai-writing-is-boring-and-dangerous/4414930), which is highly detrimental when writing tutorials.
Diverse authors bring diverse viewpoints which we may not have considered and value in our didactics.
LLMs by definition are homogeneous and homogenize all of the text they edit.

### GTN is all about community, we want you to learn and grow!

If you are a **new contributor or not a native English speaker**, we understand you may be less confident about your contribution, but our human reviewers are here to guide you through the steps to improve both your contribution and your own skills and understanding. The GTN community is very international. Most contributors are not native English speakers. It is ok for GTN content not to be perfect English. Reviewers will check for typos and grammar mistakes. We much prefer your own (possibly imperfect) words, to LLM generated text.

We fully expect that new contributors will be less confident in their contributions, and our guidance to them is to start with small contributions that they can fully understand to build confidence. We aspire to be a **welcoming community that helps new contributors grow their expertise**, but learning involves taking small steps, getting feedback, and iterating. Using an LLM [harms the learning process](http://arxiv.org/abs/2601.20245) by which novices become familiar with our framework and develop skills, without significant evidence of efficiency gains.

Passing maintainer feedback to an LLM doesn’t help anyone grow, and does not sustain our community.

### The GTN is about human interaction and collaboration

Every discussion, issue, and pull request is read and reviewed by humans.
It is a **boundary point at which people interact with each other and collaborate**.
As a community, we highly value this boundary as a place of collaboration and human connection.
Low-effort, unqualified work disrespects this boundary and puts the burden of validation solely on the maintainer.
We, as an open and inclusive community, would like to keep **our goal of growing and learning together** and from each other.

### LLMs increase the burden on our volunteer maintainers and reviewers

Contributors making use of generative AI tools oftentimes **greatly increase the maintainer burden** by not fully understanding the contents of their own contribution, not describing their contribution accurately/usefully in the PR/issue/commit message, not following up on review comments, and/or not thoroughly testing tutorials and ensuring it works as described.

### Why a technical codebase exception?

While we have strong beliefs around prose text and the deleterious effects of LLMs on prose, we have come to accept that the framework (`bin/`, `_plugins/`) is potentially different. Different to the prose, as well as to other codebases. It is rarely read by humans, it is of low quality to start with, it is not meant to be reusable or reused, it is for all intents and purposes an art project showcasing how much can be accomplished with how little, how few resources, how little investment, making the bare minimum change required to achieve a goal.

We have received a handful of LLM supported contributions to this portion of the repository which have not challenged our viewpoints significantly due to how small the patches were. Our ethical, moral, and legal issues with the use of this technology remain, and in this area we are complicit in tolerating their use.

### Are you an LLM user?

Then an LLM-free GTN benefits you!
The GTN tutorials serve as a high-quality information source for Galaxy AI solutions such as the Galaxy chatbot.
As such, to prevent model collapse, we must prevent the use of LLM generated text in the GTN.

## Societal Motivation

AI tools have numerous ethical issues. [Modern AI tools are unlawful by design](https://www.amnesty.org/en/documents/pol40/0996/2026/en/) per human rights organisation Amnesty International. [They require an unreasonable amount of energy](https://www.theverge.com/2024/5/15/24157496/microsoft-ai-carbon-footprint-greenhouse-gas-emissions-grow-climate-pledge), and [water](https://restofworld.org/2024/data-centers-environmental-issues/)  to build and operate, their models are built in part with heavily [exploited workers in unacceptable working conditions](https://time.com/6247678/openai-chatgpt-kenya-workers/), usually [without consent from or compensation for creators of the source material](https://variety.com/vip-special-reports/ai-training-consent-content-licensing-special-report-1236307878/). Furthermore, they are being used to [undermine labor and justify layoffs](https://www.forbes.com/sites/richardnieva/2025/07/17/ai-tech-layoffs/). Even if one accepts the premise they are useful, they bias [high income countries and contribute to a widening skills gap](https://doi.org/10.1126%2Fscience.adz9311). And generally are contrary to the [fundamental reasons we do science](https://terrytao.wordpress.com/2026/10/15/why-are-we-even-doing-all-this/) and we share our knowledge with others.

These are harms that we do not want to perpetuate, even if only indirectly.


## Acknowledgements

This AI Policy is heavily inspired by:

- [LLVM AI Tool Use Policy](https://llvm.org/docs/AIToolPolicy.html)
- [The Carpentries Generative AI Contributions Policy](https://docs.carpentries.org/policies/genai-policy.html)
- [rust-lang/rust LLM policy](https://blog.rust-lang.org/inside-rust/2026/08/05/rust-langrust-is-adopting-an-llm-policy/)


## Frequently Asked Questions

### Are LLMs used for review?

No; all contributions are peer reviewed. LLMs definitionally cannot be a human's peer, they would not be sufficient for peer review.

### I disagree with this policy

If you think you have a good case for why it should be changed, you may make it to us via github issue.
Human authorship is a point about which we feel very strongly.
