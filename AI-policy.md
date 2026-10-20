# GTN Generative AI Contributions Policy

The GTN is committed to providing high-quality human-written and human-curated tutorials. As such, the use of generative AI tools to create GTN content is not allowed.


## Rules

- It is **not permitted** to use generative AI to:

  - create GTN content (tutorials, faqs, images, news posts, and any other user-facing content)
  - communicate about your contribution (GitHub PR descriptions, comments, issues, reviews, chat messages, etc)

- It's fine to use LLMs to answer questions, analyze, distill, refine, check, suggest, review. But all GTN contributions and communication must be in your own words.
- For contributions to the **codebase** (`bin`, `_plugins`) you may be *assisted* by generative AI, under the following conditions:

  - You must **understand and take responsibility** for every line of code in your contribution.
  - You must be able to **answer technical questions in review** and make any requested changes.
  - You must **disclose and describe** the use of generative AI in your PR description.
  - You must write code comments and documentation yourself. Describing not only what the code does, but why you chose to do it this way.
  - You must write the PR/issue descriptions and any other (GitHub) communications yourself.

- Every **contributor is responsible** for ensuring that their contribution does not violate our policy or license
- Maintainers are encouraged to deprioritise or reject contributions that they suspect have been made in violation of this policy.


## Motivation

### The GTN is for Humans

Learners already have the ability to ask LLMs to create tutorials for them; The GTN is a resource for human written, curated and verified tutorials.

### We value diverse voices

LLMs [distort our written language](https://arxiv.org/abs/2603.18161) and result in [semantic ablation](https://www.theregister.com/software/2026/02/16/semantic-ablation-why-ai-writing-is-boring-and-dangerous/4414930), which is highly detrimental when writing tutorials.
Diverse authors bring diverse viewpoints which we may not have considered and value in our didactics.
LLMs by definition are homogeneous and homogenize all of the text they edit.

### Contribute what you can, how you can, we're happy to have you join and want you to learn and grow!

If you are a **new contributor or not a native English speaker**, we understand you may be less confident about your contribution, but our human reviewers are here to guide you through the steps to improve both your contribution and your own skills and understanding. The GTN community is very international. Most contributors are not native English speakers. It is ok for GTN content not to be perfect English. Reviewers will check for typos and grammar mistakes. We much prefer your own (possibly imperfect) words, to LLM generated text.

We fully expect that new contributors will be less confident in their contributions, and our guidance to them is to start with small contributions that they can fully understand to build confidence. We aspire to be a **welcoming community that helps new contributors grow their expertise**, but learning involves taking small steps, getting feedback, and iterating.

Passing maintainer feedback to an LLM doesn’t help anyone grow, and does not sustain our community.

### The GTN is an open, inclusive community

Every discussion, issue, and pull request is read and reviewed by humans.
It is a **boundary point at which people interact with each other and collaborate**.
As a community, we highly value this boundary as a place of collaboration and human connection.
Low-effort, unqualified work disrespects this boundary and puts the burden of validation solely on the maintainer.
We, as an open and inclusive community, would like to keep **our goal of growing and learning together** and from each other.

### LLMs have negative impacts on our work

Contributors making use of generative AI tools oftentimes **greatly increase the maintainer burden** by not fully understanding the contents of their own contribution, not describing their contribution accurately or usefully in the PR/issue/commit message, and/or not thoroughly testing tutorials and ensuring it works as described.

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
