using MPIFiles
using Base.Test
using Unitful

rG = RegularGrid{typeof(1.0u"mm")}([3,3,3],[3.0,3.0,3.0]u"mm",[0.0,0.0,0.0]u"mm")



randPos, randPosTime, numFramesToAcquire, seed = createRandPositions([10.0u"mm"; 10.0u"mm"; 10.0u"mm"], UInt64(100), UInt64(100))
randPos, randPosTime, numFramesToAcquire, seed = createRandPositions([10.0u"mm"; 10.0u"mm"; 10.0u"mm"], UInt64(100), UInt64(100),
                                                  offset =[10.0u"mm";10.0u"mm";1.0u"mm"], waitLingerTime=0.5u"s",seed=UInt64(10))

# create custom system function measurement positions
numFrames=80 # recommendation 130, if you want 100 frames, 2*15 to cut off later
coords,timeSpan = createCartesianSF([10,9,8], 10, MPILib.timeSpanPerFrame3D*numFrames)
