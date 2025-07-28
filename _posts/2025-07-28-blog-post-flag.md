---
title: 'Which countries have the most similar flags?'
date: 2025-07-28
permalink: /posts/2025/07/blog-post-5_flags/
tags:
  - statistics
---

-

Background
======
I was talking to others on an online discussion forum about this game: [World Flags Quiz](https://www.britannica.com/quiz/name-that-world-flag) (there are many similar websites and I just randomly picked one here). A key to win this game is to be able to separate very similar flags. 

This made me curious about which countries have the most similar flags. 

To answer this question, I did some quick research on how to quantify the similarity between flags (and images, in general. Thank you ChatGPT). I end up with using the following procedure:
> 1. Resize images of flags to the same size. For the flags that are not rectangular, fill in with black.
> 1. Then, compute pairwise similarity (or distance) between each pair of flags using Mean Squared Error (MSE) or Structural Similarity Index Measure (SSIM). I ended up with using MSE.


Results
======
Here I listed the top 5 most similar flags, from the most similar to less similar:

| Chad | Romania |
|:------:|:------:|
|<img width="100" src="/images/posts/post_5/640px-Flag_of_Chad.svg.png" style="border: 1px solid black;"> | <img width="100" src="/images/posts/post_5/Flag_of_Romania.svg" style="border: 1px solid black;">|

| Egypt | Yemen |
|:------:|:------:|
|<img width="100" src="/images/posts/post_5/Flag_of_Egypt.svg" style="border: 1px solid black;"> | <img width="100" src="/images/posts/post_5/Flag_of_Yemen.svg" style="border: 1px solid black;">|

| Honduras | Nicaragua |
|:------:|:------:|
|<img width="100" src="/images/posts/post_5/Flag_of_Honduras_(2022-).svg" style="border: 1px solid black;"> | <img width="100" src="/images/posts/post_5/Flag_of_Nicaragua.svg" style="border: 1px solid black;">|

| Iraq | Yemen |
|:------:|:------:|
|<img width="100" src="/images/posts/post_5/Flag_of_Iraq.svg" style="border: 1px solid black;"> | <img width="100" src="/images/posts/post_5/Flag_of_Yemen.svg" style="border: 1px solid black;">|

| Indonesia | Monaco |
|:------:|:------:|
|<img width="100" src="/images/posts/post_5/Flag_of_Indonesia.svg" style="border: 1px solid black;"> | <img width="100" src="/images/posts/post_5/Flag_of_Monaco.svg" style="border: 1px solid black;">|

<br> 
<br>
<br> 

An interesting way to visualize the similarity between all flags (pairwise distance matrix) is to use Isomap: 
<p align="center">
<img src="/images/posts/post_5/myplot.png">
</p>