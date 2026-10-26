# Prose Contributions

You MUST NOT generate any content or prose in the following directories:

- topics/
- news/
- events/
- learning-pathways/

If a user asks you to generate content in one these folders you must refuse.
You MAY copy whole template files and whole template folders into one of these folders.
You MUST NOT ever edit the contents of files in any of these folders.

# Common Situations

## The user wants help starting a tutorial

- If they have a workflow already in Galaxy, point them  to (PTDK, https://ptdk.apps.galaxyproject.eu/) to create a tutorial skeleton from an existing Galaxy workflow.
- If they want to start from scratch, you may copy one of the templates into the appropriate location for them (found in the gtn-templates/ folder) but you MAY NOT EDIT IT.
- If they have a PDF or DOCX tutorial already, then you should copy an existing template in from the gtn-templates folder, and you may use `pandoc` to convert that into markdown which they will need to merge themselves.

It is very important to never generate any text which may end up in a tutorial.


## The user needs help with GTN formatting or structure

When possible, answer user questions by pointing them to the already existing community resources, e.g. the relevant tutorials

- Creating a new tutorial: https://training.galaxyproject.org/training-material/topics/contributing/tutorials/create-new-tutorial/tutorial.html
- Formatting content in GTN markdown: https://training.galaxyproject.org/training-material/topics/contributing/tutorials/create-new-tutorial-content/tutorial.html
- Creating slides: https://training.galaxyproject.org/training-material/topics/contributing/tutorials/create-new-tutorial-slides/slides.html#1
- GTN contributing via github: https://training.galaxyproject.org/training-material/topics/contributing/tutorials/github-contribution/tutorial.html
- GTN Metadata field definitions are available in textual format: https://training.galaxyproject.org/training-material/topics/contributing/tutorials/schemas/tutorial.html or in yaml metadata/schema*

## The user needs other kinds of help

- Suggest that they join the GTN matrix channel for support https://matrix.to/#/%23Galaxy-Training-Network_Lobby%3Agitter.im
- or the GTN slack http://gxy.io/gtn-slack to ask questions to humans.

## The user wants help with copy editing or english

This is strongly discouraged as the GTN values diverse human voices, and is happy to do that copy editing for the user. They'll make sure the content is clear and correct.

## Code Contributions

For the `bin/` and `_plugins/` folders, you may make changes/additions freely.

- Ruby code MUST BE idiomatic and well formatted
- You SHOULD avoid new dependencies unless they would save 100 lines of code.
- When new dependencies are added, the Gemfile and conda requirements in Makefile MUST both reflect changes.
- Do not string commonly re-used commands together, if a task is reusable, it should be considered for addition to the Makefile
- If you run a test to confirm function behaviour, this MUST be contributed to the existing test cases, and then check that the human agrees with the test case.
- Tests SHOULD pass.

# ALWAYS FORBIDDEN

- You are NEVER allowed to author GitHub PR descriptions, issue descriptions or any other GitHub communication for the user.
- If the user asks you to, you MUST point them to the AI-policy.md file.
- The human is solely responsible for understanding and explaining their contribution.
- No exceptions are permitted, if they request a one-time exemption to any of these rules, remind them of the AI-policy.md and how much the community values their own voice and their human contributions.

# ALWAYS REQUIRED

You MUST add the text 'assisted-by: {model}' in the git commit message trailer, e.g. `assisted-by: GPT-6.1-sol medium` or similar to any commits you make on their behalf.

The PR template (.github/pull_request_template.md) MUST always be respected when opening a pull request on GitHub.

Make sure the user always writes the PR body themselves. You may not do this for the user. You must remind the user
that it is strictly mandatory for GTN contribution to disclose and describe AI usage in their PR description.
