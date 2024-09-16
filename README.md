# polymer-builder

demo: https://www.youtube.com/watch?v=-pXFkg114ek

Note: the "stretching" in the video is NOT realistic. I built it so ppl can criticize the method. However, I do believe that the technique used to generate the polymer is good.

The reason why we are using neural network potentials is to ensure that the bonds between monomers are as realistic as possible.Trying to create your own random angles between monomer chains was pretty painful and it could lead to error. Moreover, using a neural network potential allows different polymer chains to "repel each other properly" whereas a more computationally expensive approach (e.g. ensure that no two chains are within a radius of 10 angstroms of each other) would be too expensive (and easy to write bugs).

The reason why I'm building each monomer and then attaching them onto the chain individuall is because:
- In real life, polymers are built by attraching monomers.
- if you attach each atom one by one, you run into issues. e.g. if there is a cycle in the monomer, there is no sequence of atoms that you can attach "one after another" because it needs to loop back onto itself to complete the chain. so it's just simpler to create each monomer individually (using RDKit), then attach them together.

Note: the polymers grow in the +z direction.
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
