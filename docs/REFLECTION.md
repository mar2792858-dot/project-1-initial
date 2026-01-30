# Reflection Log

This document captures reflections on the development of 3D geometric classes in Java, focusing on design patterns, principles, and lessons learned.

I researched Composite design pattern and utilized this website:
https://www.geeksforgeeks.org/java/composite-design-pattern-in-java/.

This pattern seems very similar to the parent-child inheritance in C++
where the child class utilizes methods from the parent class, but
can override or have its own. The original class methods do not change,
but the objects that are created have the ability to set and change
characteristics about the object. I think this pattern is seen in the
Cube3D class where the constructor has the access level (public) the
constructor name, and the argument passed in is the vertices from the 
Point3D class. It also uses the distanceTo function of the Point3D class.
There is validation of the vertices and a copy is made as well. The
constructor ends with a logger. This pattern is seen in the other methods
of the Cube3D class where there is a characteristic from another class
either being passed into the argument or used within that method. 

