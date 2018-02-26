Here are some brief notes explaining the pictures in this file.

There are 3 different datasets used here. They are called DFC (the urban scene), Umbrella, and Jadeplant.
For each dataset we use the weighting from the 2017 ICIP paper, along with the Nystrom method to get 100 eigenvectors.
This section of the code was done in Matlab. Those 100 eigenvectors are used both in the semisupervised MBO algorithm,
as well as the spectral clustering algorithm.

This file contains the 100 eigenvectors associated with each dataset (one folder per set of eigenvectors)
as well as the input to and results from each algorithm (one folder per set of results)


Umbrella:
Here we segment into 6 classes.
1. front umbrella
2. back umbrella
3. brown cabinet
4. white background wall
5. gray background wall
6. miscellaneous darker background object
The fidelity consists of small rectangles highlighting important objects, and consists of about 5% of the
total image. Overall, this worked well with very little tinkering around


DFC:
Here we segment into 6 classes. 
1. Highest roofs
2. Mid-level roofs
3. Light gray road (in the upper right)
4. Dark gray road (between the houses)
5. Grass
6. Bright white objects
The algorithm had a very hard time differentiating the grass from the dark gray road. Even in the final results
you can see some misclassification here. To handle this I put more effort into choosing the fidelity (not just rectangles).
The grass was chosen based on color, and the roads were represented by lines. The first fidelity input that gave a good result
used 18% of the total pixels. From this I got down to 5% by randomly removing 13/18 of the fidelity pixels. We should probably
talk about whether or not this is reasonable. To me it feels like "gaming" the results.


Jadeplant:
Here we segment into 9 classes.
1. Dark black background
2. Brown box (near the back)
3. Wooden table (in the foreground. Everything is resting on this table)
4. Red pot (that the plant is in)
5. Yellow cloth (behind the plant)
6. Plant branches
7. Plant leaves
8. The blue grid thing at the front
9. The small wooden block that the blue grid rests on
Once again, the fidelity consists of 5% of the pixels. I chose easy rectangles whenever they were available.
For the leaves I chose points based on color. And same for the grid. The plant branches were the hardest to pick a fidelity for,
and I think the end result shows this (branches, leaves, and the blue grid are all sort of merged). I'm sure I could make it work
by once again very carefully choosing the fidelity points, but like the DFC case, this seems disingenuous.
Here, the algorithm had a very difficult time separating the blue grid from the red pot (which I don't understand,
since the colors are very different). This is why the fidelity contains pretty much the entire blue grid.
