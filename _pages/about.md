---
layout: about
title: about
permalink: /
subtitle: Mathematical modeler

profile:
  align: right
  image: covid_ts.png
  image_circular: false # crops the image to make it circular
  address: >
    <p>Disease incidence and contact networks</p>
    <p> </p>

news: false  # includes a list of news items
selected_papers: true # includes a list of papers marked as "selected={true}"
social: true  # includes social icons at the bottom of the page
announcements:
  enabled: false # includes a list of news items
  scrollable: true # adds a vertical scroll bar if there are more than 3 news items
  limit: 5 # leave blank to include all the news in the `_news` folder
latest_posts:
  enabled: false
  scrollable: true # adds a vertical scroll bar if there are more than 3 new posts items
  limit: 3 # leave blank to include all the blog posts
---

{% assign currentYear = "now" | date: "%Y" | plus: 0 %}
{% assign startYear = 2012 %}
{% assign diff = currentYear | minus: startYear %}

A mathematician by [trade](http://tsp.imath.kiev.ua/files/242/tsp1710_12.pdf), I have been working in **disease modeling**  for the past {{ diff }} years. 


Infectious disease dynamics are **multi-scale** by definition. The combination of **pathogen** or **immune dynamics** within a host and changing patterns of **interactions between hosts** during transmission result in rich population-level phenomena. Ranging from stochastic emergence and extinction, to structured and well-mixed epidemic processes, I am interested in applying **mathematical** and **statistical techniques** to answer questions from the angle of **population health**:  

* What is the role of host response to infection towards [disease outcome](https://doi.org/10.1186/1471-2164-15-1161)?

* When does within-host [heterogeneity](https://doi.org/10.1098/rspb.2022.0232) matter for disease dynamics at the population level?

* How can we [disentangle](https://doi.org/10.1016/j.ijid.2021.08.029) various [sources](https://doi.org/10.1371/journal.pbio.2002468) of heterogeneity to evaluate [intervention](https://doi.org/10.1016/j.vaccine.2019.02.073) effectiveness?

I have been mainly interested in pathogens and diseases affecting humans such as HIV, Influenza, Polio, HPV, Coronaviruses and Malaria.

I am currently employed at the [Swiss Tropical and Public Health Institute](https://swisstph.ch).

For more details, you can download my [resume](/assets/pdf/resume_cselinger_20260324.pdf).
