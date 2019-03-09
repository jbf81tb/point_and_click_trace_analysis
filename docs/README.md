# Analyze clathrin-coated pits with mouse and keyboard input
Users are presented with a large window containing an image of the cell. The user can scroll along the frame slider to search through the movie for an ideal structure. If the user considers the features too small to see very well, they can zoom in on a quadrant for a 4x apparent increase in size. Once an interesting structure is located the user should adjust the slider near the beginning of the lifetime of the structure and then click on it.

Upon clicking, the rest of the windows fill with information on the currently selected structure. Two windows are devoted to increasing zooms of the structure selected. Every window will automatically center on the structure. 

On the most zoomed image will be circles corresponding to the automatically calculated area of the structure. Occasionally, the structure will have a "C" shape and so the area calculation will not be able to close the circle. In this case the user can left-click on this window to select a pixel that should be considered in the area, and once the user has closed the circle the area will automatically fill in. For dim structures the algorithm might consider too much to be in the area so a right-click can deselect unnecessary pixels.

The second-most zoomed image allows the user to have a wider view of neighboring structures while keeping a tight view on the structure of interest. In this way they can make sure that no other structures interfere with the structure of interest.

Additionally, the user is provided with a graph of the intensity from the TIRF image of the cell. The real data is presented with points and the best fit of a 2D Gaussian is displayed as a surface.

Once the object is selected, the user can navigate the frames using the keyboard, with the goal of locating the very first frame of the structure's existence. A keyboard command is used to bookmark this frame. Next, the user should used the keyboard to navigate to the final frame of the structure's existence and bookmark that frame. Once a beginning and end have been bookmarked, the user can use a keyboard command to save the trace into the trace database. Once saved, the structure will be circled on the main image of the cell so it doesn't accidentally get visited again.

As the user proceeds through the frames, windows in the lower right corner will track the intensity and area of frames visited so far. The user can also click on a trace that has already been saved and that will fill in these graphs with information from that cell.

When exploring a trace that has already been saved the window displaying the intensity will be disabled, because the 2D-Gaussian fit is computationally intensive and slows down the user's ability to quickly scrub through the movie.
