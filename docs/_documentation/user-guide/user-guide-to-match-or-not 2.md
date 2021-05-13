---
layout: default
title: Whether to use Matching
nav_order: 2
permalink: /user-guide/to-match-or-not
parent: User Guide
---

# To Match or Not To Match
{: .no_toc }
That is the question

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>


## Determining Whether to Use Matching Methods

Matching of treatment and control units can be a good method in order to determine treatment effects. However, certain criteria must be upheld in order for matching to be an appropriate solution for a given dataset. If these criteria are not upheld, perhaps other approaches to causal inference should be used in place of, or in addition to matching. 

### The Stable Unit Treatment Value Assumption (SUTVA)

Treatments applied to one unit should not affect the outcome of another unit. Units can not interfere with one another. This is reasonable in many situations: If two individuals are not in contact with each other, how would one individual taking a pain medication impact the outcome of another individual. 

We should also assume that the treatment doesn't have varying forms, and is completely binary. Individuals can not have taken pain medication of different strengths. 

### The Unconfoundedness Assumption

This is also referred to as "ignorability". It is important that the outcome is independent of the treatment when observable covaraiates are held constant. Omitted variable bias is a common issue that occurs when a variable impacts both treatment and outcomes, and appears in a bias of treatment effect estimates. 

In the example about pain medications, if a researcher fails to include in their dataset some underlying health condition that impacts response to pain medication, the impact of taking pain medication for a headache might be evaluated incorrectly.


### Overlap of Treatment and Control Groups

A common problem in causal inference is overlap or imbalance between treatment and control groups. A treatment and control group would have no overlap if none of the covariates have the same values. In this case, the FLAME and DAME algorithms would not find any matches, and no treatment effect estimates would be possible. 


<div class="language-markdown highlighter-rouge">
  <h4>Further Readings</h4>
  <br/>
  For more information on causal inference research and its assumptions and issues, we recommend 
  <a href="https://books.google.com/books?hl=en&lr=&id=Bf1tBwAAQBAJ&oi=fnd&pg=PR17&ots=jeVGafZSDE&sig=x9LYF4V9-wYQRQRxpudyA-d9nI0">
    Imbens, Guido W., and Donald B. Rubin. <i>Causal inference in statistics, social, and biomedical sciences</i>.
  </a>
</div>
