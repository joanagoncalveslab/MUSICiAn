# Ideas

This document is to keep track of different ideas throughout the project so that nothing gets lots


## Models

#### Attention mechanism-based model

Right now, we are focused on finding which categories lead to the best predictions. Perhaps we can use a deep learning attention mechanism based method to remove the guesswork and the need to encode repair outcome profiles. We just provide all the visible repair outcomes at once:
i.e an n x m matrix, where n is the repair outcomes and m is a vector describing the repair outcome (deletion, length, start location)