# An adaptive feature-preserving isotropic remesher

![alt text](https://github.com/rozentill/Remesher/raw/master/fig/adaptive.png "")


This is an implementation of modified classical isotropic remeshing algorithm in libigl.

#### How to build:
- Use the files under ```src```
- Set 3rd party libraries(libigl and the ```external``` in it)

#### How to run:
- A prebuilt demo is provided in ```exe```, run it by specify the model:
> remesher.exe ../data/cow_head.obj

#### How to use:
- 'Reset' firstly.
- Set target length, mostly 0.02 (for cactus), 0.01 (for horse) or 0.005 (for cow_head)
- Click 'Remesh' to do one iteration
- Do at least 4 iterations
- Similary, you can click 'adaptive', now only tried on cow_head (0.008 length)
