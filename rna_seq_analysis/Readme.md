# Overview

Here are the RNA-Seq analysis pipeline templates that I used.

The main pipeline that I utilize is the New-Tuxedo suite 'hisat2 + stringtie + ballgown'. I am genuiely surprised by how fast tihs pipeline works compared to the speed of the old suite 'bowtie + tophat + cufflinks'.

Again, I wish that my project could move forward to the Drop-Seq so that I can tease out the molecular background of the heterogeneous iEC/ iEC daughter cells. To prepare for that, I played with RNA-Seq in 92 yeast strains.

Another part is using R with Juputer notebook for RNA-Seq data analysis. Here is a quote about the benefit of doing this:

> Jupyter was designed to enable sharing of notebooks with other people. The idea is that you can write some code, mix some text with the code, and publish this as a notebook.  In the notebook they can see the code as well as the actual results of running the code.

> This is a nice way of sharing little experimental snippets, but also to publish more detailed reports with explanations and full code sets.  Of course, a variety of web services allows you to post just code snippets (e.g. [gist](https://gist.github.com/)). What makes Jupyter different is that the service will actually render the code output.

> One interesting benefit of using Jupyter is that Github magically renders notebooks.