# polymer-builder

demo: https://www.youtube.com/watch?v=-pXFkg114ek

Note: the "stretching" in the video is NOT realistic. I built it so ppl can criticize the method. However, I do believe that the technique used to generate the polymer is good.

The reason why we are using neural network potentials is to ensure that the bonds between monomers are as realistic as possible.Trying to create your own random angles between monomer chains was pretty painful and it could lead to error. Moreover, using a neural network potential allows different polymer chains to "repel each other properly" whereas a more computationally expensive approach (e.g. ensure that no two chains are within a radius of 10 angstroms of each other) would be too expensive (and easy to write bugs).


In the future, we can have realtime MD so we can drag these polymers around like a ragdoll and attach cross link bonds between polymer chains.

### how to use the viewer

go into viewer, then run these commands:
```
npm install
npm run dev
```
then go into http://localhost:3000/ on your browser!


### How to generate the polymer:
`python assemble.py`

### How to stretch the polymer (not realistic):
`python stretch.py`
