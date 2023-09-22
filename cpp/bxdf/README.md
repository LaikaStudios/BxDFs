# BxDF Plugin Notes

While generating the [`doc/media`](../../doc/media) images at HD resolution,
I found that [`PraterFuzz`](./PraterFuzz.cpp) causes PRMan to crash 
upon render completion.
This doesn't happen with smaller resolution renderings or when Live Rendering.
My guess is it has something to do with the fairly large globally declared response sampling direction vector array,
but I have no idea why that would cause a problem.
If anyone can provide some insight that would be appreciated.

As mentioned in the docs, I do plan on eliminating the use of this array and replacing the numerical normalization that utilizes it with an approximate analytical normalization function.
I'm currently waiting for LAIKA to purchase a Matlab license so I can proceed with that effort.
