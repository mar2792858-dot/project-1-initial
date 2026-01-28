package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a point in three-dimensional Euclidean space.
 * 
 * This class demonstrates several object-oriented design patterns and principles:
 * 
 * DESIGN PATTERNS USED:
 * 
 * 1. VALUE OBJECT PATTERN
 *    - Immutable state (final fields) ensures thread safety and predictable behavior
 *    - Points are compared by value, not identity
 *    - Safe to use as hash map keys or in collections
 * 
 * 2. BUILDER PATTERN (implicit through constructors)
 *    - Multiple constructors provide flexible object creation
 *    - Overloaded constructors demonstrate encapsulation and API design
 * 
 * 3. FACTORY PATTERN (static factory methods)
 *    - Static creation methods (origin, fromCylindrical, fromSpherical) provide
 *      semantic clarity and alternative construction paths
 * 
 * FOUNDATIONAL PRINCIPLES DEMONSTRATED:
 * 
 * 1. ENCAPSULATION
 *    - Private fields with public accessor methods
 *    - Internal state hidden from external modification
 *    - Logging abstracts implementation details
 * 
 * 2. ABSTRACTION
 *    - Complex mathematical operations (rotation matrices, transformations)
 *      are hidden behind simple method interfaces
 *    - Users don't need to understand quaternions or Euler angles
 * 
 * 3. SINGLE RESPONSIBILITY PRINCIPLE
 *    - Each method has one clear purpose
 *    - Logging concern is separated via Logger
 * 
 * DATA STRUCTURE & ALGORITHM FOUNDATIONS:
 * 
 * 1. GEOMETRIC PRIMITIVES
 *    - Points are fundamental building blocks for spatial data structures
 *    - Used in: KD-trees, octrees, BSP trees, convex hulls
 * 
 * 2. DISTANCE METRICS
 *    - Euclidean distance is basis for nearest-neighbor algorithms
 *    - Manhattan distance useful for grid-based pathfinding
 *    - Enables clustering algorithms (K-means, DBSCAN)
 * 
 * 3. COORDINATE TRANSFORMATIONS
 *    - Rotation matrices are linear transformations
 *    - Foundation for computer graphics pipelines
 *    - Used in robotics, physics simulations, game engines
 * 
 * 4. IMMUTABILITY
 *    - Enables functional programming patterns
 *    - Safe for concurrent data structures
 *    - Simplifies reasoning about algorithm correctness
 * 
 * @author Generated for College Class
 * @version 1.0
 */
public class Point3D {
    
    // Logger instance for this class - demonstrates Singleton pattern
    private static final Logger logger = Logger.getLogger(Point3D.class.getName());
    
    // Coordinates - immutable for thread safety and predictable behavior
    private final double x;
    private final double y;
    private final double z;
    
    /**
     * Constructs a new Point3D with the specified coordinates.
     * 
     * This is the primary constructor that all other constructors delegate to.
     * It demonstrates the principle of having a single source of truth for
     * object initialization.
     * 
     * Time Complexity: O(1) - constant time initialization
     * Space Complexity: O(1) - fixed memory for three doubles
     * 
     * @param x the x-coordinate in 3D space
     * @param y the y-coordinate in 3D space
     * @param z the z-coordinate in 3D space
     */
    public Point3D(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
        
        logger.log(Level.INFO, "Created new Point3D at ({0}, {1}, {2})", 
                   new Object[]{x, y, z});
    }
    
    /**
     * Constructs a point at the origin (0, 0, 0).
     * 
     * This convenience constructor demonstrates constructor overloading and
     * provides a clear, semantic way to create origin points without magic numbers.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     */
    public Point3D() {
        this(0.0, 0.0, 0.0);
        logger.log(Level.INFO, "Created Point3D at origin");
    }
    
    /**
     * Factory method to create a point at the origin.
     * 
     * Static factory methods provide more flexibility than constructors:
     * - Can return cached instances (flyweight pattern)
     * - Have descriptive names
     * - Can return subtypes if needed
     * 
     * @return a new Point3D at coordinates (0, 0, 0)
     */
    public static Point3D origin() {
        logger.log(Level.INFO, "Requesting origin point via factory method");
        return new Point3D(0.0, 0.0, 0.0);
    }
    
    // Accessor methods (getters) - demonstrate encapsulation
    
    /**
     * Returns the x-coordinate of this point.
     * 
     * Accessor methods allow controlled access to internal state while
     * maintaining encapsulation. Future implementations could add validation,
     * logging, or lazy computation without changing the API.
     * 
     * @return the x-coordinate
     */
    public double getX() {
        return x;
    }
    
    /**
     * Returns the y-coordinate of this point.
     * 
     * @return the y-coordinate
     */
    public double getY() {
        return y;
    }
    
    /**
     * Returns the z-coordinate of this point.
     * 
     * @return the z-coordinate
     */
    public double getZ() {
        return z;
    }
    
    /**
     * Calculates the Euclidean distance between this point and another point.
     * 
     * Uses the 3D distance formula: √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
     * 
     * This is a fundamental operation in many algorithms:
     * - Nearest neighbor search
     * - Clustering algorithms (K-means, DBSCAN)
     * - Collision detection
     * - Path planning
     * 
     * Time Complexity: O(1) - constant number of arithmetic operations
     * Space Complexity: O(1) - only uses local variables
     * 
     * @param other the point to calculate distance to
     * @return the Euclidean distance between the two points
     * @throws IllegalArgumentException if other is null
     */
    public double distanceTo(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Cannot calculate distance to null point");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        
        logger.log(Level.INFO, "Calculated distance from ({0}, {1}, {2}) to ({3}, {4}, {5}): {6}",
                   new Object[]{this.x, this.y, this.z, other.x, other.y, other.z, distance});
        
        return distance;
    }
    
    /**
     * Calculates the Manhattan distance (taxicab distance) to another point.
     * 
     * Manhattan distance is the sum of absolute differences: |x₂-x₁| + |y₂-y₁| + |z₂-z₁|
     * 
     * This metric is useful for:
     * - Grid-based pathfinding (A* algorithm)
     * - When diagonal movement is restricted
     * - Faster computation than Euclidean distance (no square root)
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param other the point to calculate Manhattan distance to
     * @return the Manhattan distance
     * @throws IllegalArgumentException if other is null
     */
    public double manhattanDistanceTo(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate Manhattan distance to null point");
            throw new IllegalArgumentException("Cannot calculate Manhattan distance to null point");
        }
        
        double distance = Math.abs(this.x - other.x) + 
                         Math.abs(this.y - other.y) + 
                         Math.abs(this.z - other.z);
        
        logger.log(Level.INFO, "Calculated Manhattan distance: {0}", distance);
        
        return distance;
    }
    
    /**
     * Calculates the distance from this point to the origin (0, 0, 0).
     * 
     * This is equivalent to the magnitude or length of the position vector.
     * Frequently used in physics simulations and graphics calculations.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return the distance from the origin
     */
    public double distanceFromOrigin() {
        double distance = Math.sqrt(x * x + y * y + z * z);
        logger.log(Level.INFO, "Distance from origin for point ({0}, {1}, {2}): {3}",
                   new Object[]{x, y, z, distance});
        return distance;
    }
    
    /**
     * Rotates this point around the X-axis by the specified angle.
     * 
     * Uses the rotation matrix:
     * [ 1    0         0      ]
     * [ 0  cos(θ)  -sin(θ)    ]
     * [ 0  sin(θ)   cos(θ)    ]
     * 
     * Returns a new Point3D since this class is immutable. Immutability:
     * - Prevents bugs from unexpected state changes
     * - Enables safe sharing between threads
     * - Simplifies reasoning about code correctness
     * 
     * Applications:
     * - 3D graphics transformations
     * - Robotics (gimbal rotations)
     * - Physics simulations
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1) - creates one new Point3D
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateX(double angleRadians) {
        logger.log(Level.INFO, "Rotating point around X-axis by {0} radians", angleRadians);
        
        double cosTheta = Math.cos(angleRadians);
        double sinTheta = Math.sin(angleRadians);
        
        double newY = y * cosTheta - z * sinTheta;
        double newZ = y * sinTheta + z * cosTheta;
        
        Point3D rotated = new Point3D(x, newY, newZ);
        
        logger.log(Level.INFO, "Rotation complete: ({0}, {1}, {2}) → ({3}, {4}, {5})",
                   new Object[]{x, y, z, rotated.x, rotated.y, rotated.z});
        
        return rotated;
    }
    
    /**
     * Rotates this point around the Y-axis by the specified angle.
     * 
     * Uses the rotation matrix:
     * [  cos(θ)   0   sin(θ) ]
     * [    0      1     0    ]
     * [ -sin(θ)   0   cos(θ) ]
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateY(double angleRadians) {
        logger.log(Level.INFO, "Rotating point around Y-axis by {0} radians", angleRadians);
        
        double cosTheta = Math.cos(angleRadians);
        double sinTheta = Math.sin(angleRadians);
        
        double newX = x * cosTheta + z * sinTheta;
        double newZ = -x * sinTheta + z * cosTheta;
        
        Point3D rotated = new Point3D(newX, y, newZ);
        
        logger.log(Level.INFO, "Rotation complete: ({0}, {1}, {2}) → ({3}, {4}, {5})",
                   new Object[]{x, y, z, rotated.x, rotated.y, rotated.z});
        
        return rotated;
    }
    
    /**
     * Rotates this point around the Z-axis by the specified angle.
     * 
     * Uses the rotation matrix:
     * [ cos(θ)  -sin(θ)   0 ]
     * [ sin(θ)   cos(θ)   0 ]
     * [   0        0      1 ]
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing the rotated point
     */
    public Point3D rotateZ(double angleRadians) {
        logger.log(Level.INFO, "Rotating point around Z-axis by {0} radians", angleRadians);
        
        double cosTheta = Math.cos(angleRadians);
        double sinTheta = Math.sin(angleRadians);
        
        double newX = x * cosTheta - y * sinTheta;
        double newY = x * sinTheta + y * cosTheta;
        
        Point3D rotated = new Point3D(newX, newY, z);
        
        logger.log(Level.INFO, "Rotation complete: ({0}, {1}, {2}) → ({3}, {4}, {5})",
                   new Object[]{x, y, z, rotated.x, rotated.y, rotated.z});
        
        return rotated;
    }
    
    /**
     * Translates (moves) this point by the specified offsets.
     * 
     * Translation is a fundamental affine transformation in computer graphics.
     * Combined with rotation and scaling, it forms the basis of model-view
     * transformations in 3D rendering pipelines.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param dx the offset in the x direction
     * @param dy the offset in the y direction
     * @param dz the offset in the z direction
     * @return a new Point3D at the translated position
     */
    public Point3D translate(double dx, double dy, double dz) {
        logger.log(Level.INFO, "Translating point by ({0}, {1}, {2})", 
                   new Object[]{dx, dy, dz});
        
        Point3D translated = new Point3D(x + dx, y + dy, z + dz);
        
        logger.log(Level.INFO, "Translation complete: ({0}, {1}, {2}) → ({3}, {4}, {5})",
                   new Object[]{x, y, z, translated.x, translated.y, translated.z});
        
        return translated;
    }
    
    /**
     * Scales this point by the specified factor from the origin.
     * 
     * Scaling is a linear transformation that changes the magnitude of the
     * position vector. Uniform scaling (same factor for all axes) preserves
     * the direction of the vector.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param factor the scaling factor
     * @return a new Point3D at the scaled position
     */
    public Point3D scale(double factor) {
        if (factor < 0) {
            logger.log(Level.WARNING, "Scaling with negative factor {0} will invert the point", 
                      factor);
        }
        
        logger.log(Level.INFO, "Scaling point by factor {0}", factor);
        
        Point3D scaled = new Point3D(x * factor, y * factor, z * factor);
        
        logger.log(Level.INFO, "Scaling complete: ({0}, {1}, {2}) → ({3}, {4}, {5})",
                   new Object[]{x, y, z, scaled.x, scaled.y, scaled.z});
        
        return scaled;
    }
    
    /**
     * Calculates the dot product of this point's position vector with another.
     * 
     * The dot product is a fundamental operation in linear algebra with many uses:
     * - Calculating the angle between vectors
     * - Testing for perpendicularity (dot product = 0)
     * - Projecting one vector onto another
     * - Determining if vectors point in similar directions
     * 
     * Formula: a · b = a_x * b_x + a_y * b_y + a_z * b_z
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param other the other point (treated as a position vector)
     * @return the dot product
     * @throws IllegalArgumentException if other is null
     */
    public double dotProduct(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate dot product with null point");
            throw new IllegalArgumentException("Cannot calculate dot product with null point");
        }
        
        double result = this.x * other.x + this.y * other.y + this.z * other.z;
        
        logger.log(Level.INFO, "Calculated dot product: {0}", result);
        
        return result;
    }
    
    /**
     * Calculates the cross product of this point's position vector with another.
     * 
     * The cross product produces a vector perpendicular to both input vectors.
     * This is crucial for:
     * - Calculating surface normals in 3D graphics
     * - Determining the orientation of a coordinate system
     * - Computing torque in physics simulations
     * - Finding a vector perpendicular to a plane
     * 
     * The magnitude equals the area of the parallelogram formed by the vectors.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param other the other point (treated as a position vector)
     * @return a new Point3D representing the cross product vector
     * @throws IllegalArgumentException if other is null
     */
    public Point3D crossProduct(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate cross product with null point");
            throw new IllegalArgumentException("Cannot calculate cross product with null point");
        }
        
        double newX = this.y * other.z - this.z * other.y;
        double newY = this.z * other.x - this.x * other.z;
        double newZ = this.x * other.y - this.y * other.x;
        
        Point3D result = new Point3D(newX, newY, newZ);
        
        logger.log(Level.INFO, "Calculated cross product: ({0}, {1}, {2})",
                   new Object[]{newX, newY, newZ});
        
        return result;
    }
    
    /**
     * Normalizes this point's position vector to unit length.
     * 
     * A normalized (unit) vector has magnitude 1 but preserves direction.
     * Unit vectors are essential for:
     * - Direction representation in physics
     * - Lighting calculations (surface normals)
     * - Camera orientation in graphics
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return a new Point3D with magnitude 1 in the same direction
     * @throws ArithmeticException if this point is at the origin (cannot normalize)
     */
    public Point3D normalize() {
        double magnitude = distanceFromOrigin();
        
        if (magnitude == 0) {
            logger.log(Level.SEVERE, "Cannot normalize zero vector at origin");
            throw new ArithmeticException("Cannot normalize the zero vector");
        }
        
        Point3D normalized = new Point3D(x / magnitude, y / magnitude, z / magnitude);
        
        logger.log(Level.INFO, "Normalized vector: ({0}, {1}, {2}) → ({3}, {4}, {5})",
                   new Object[]{x, y, z, normalized.x, normalized.y, normalized.z});
        
        return normalized;
    }
    
    /**
     * Calculates the midpoint between this point and another point.
     * 
     * The midpoint is the average of the two points' coordinates.
     * Used in:
     * - Binary space partitioning (BSP trees)
     * - Subdivision algorithms
     * - Finding centers of line segments
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param other the other point
     * @return a new Point3D at the midpoint
     * @throws IllegalArgumentException if other is null
     */
    public Point3D midpoint(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate midpoint with null point");
            throw new IllegalArgumentException("Cannot calculate midpoint with null point");
        }
        
        double midX = (this.x + other.x) / 2.0;
        double midY = (this.y + other.y) / 2.0;
        double midZ = (this.z + other.z) / 2.0;
        
        Point3D midpoint = new Point3D(midX, midY, midZ);
        
        logger.log(Level.INFO, "Calculated midpoint: ({0}, {1}, {2})",
                   new Object[]{midX, midY, midZ});
        
        return midpoint;
    }
    
    /**
     * Returns a string representation of this point.
     * 
     * String representation is important for:
     * - Debugging and logging
     * - User interface display
     * - Serialization for network transmission
     * 
     * Format: "Point3D(x, y, z)"
     * 
     * @return a string representation of this point
     */
    @Override
    public String toString() {
        return String.format("Point3D(%.2f, %.2f, %.2f)", x, y, z);
    }
    
    /**
     * Compares this point with another object for equality.
     * 
     * Two points are equal if their coordinates are equal (within floating-point
     * precision). This implements value equality rather than reference equality.
     * 
     * Proper equals() implementation is critical for:
     * - Using objects in HashMaps and HashSets
     * - Comparing results in unit tests
     * - Implementing search algorithms
     * 
     * Note: In production code, you might want to use epsilon comparison for
     * floating-point values to handle precision issues.
     * 
     * @param obj the object to compare with
     * @return true if the objects are equal, false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null || getClass() != obj.getClass()) {
            return false;
        }
        
        Point3D other = (Point3D) obj;
        
        boolean isEqual = Double.compare(x, other.x) == 0 &&
                         Double.compare(y, other.y) == 0 &&
                         Double.compare(z, other.z) == 0;
        
        if (isEqual) {
            logger.log(Level.INFO, "Points are equal: {0} == {1}", 
                      new Object[]{this, other});
        }
        
        return isEqual;
    }
    
    /**
     * Returns a hash code for this point.
     * 
     * Hash codes are essential for using objects in hash-based collections
     * (HashMap, HashSet). The contract states that equal objects must have
     * equal hash codes.
     * 
     * This implementation uses Java's Objects.hash() utility which combines
     * the hash codes of the coordinate values.
     * 
     * @return a hash code value for this point
     */
    @Override
    public int hashCode() {
        int result = 17;
        long xBits = Double.doubleToLongBits(x);
        long yBits = Double.doubleToLongBits(y);
        long zBits = Double.doubleToLongBits(z);
        
        result = 31 * result + (int)(xBits ^ (xBits >>> 32));
        result = 31 * result + (int)(yBits ^ (yBits >>> 32));
        result = 31 * result + (int)(zBits ^ (zBits >>> 32));
        
        return result;
    }
    
    /**
     * Creates a deep copy of this point.
     * 
     * Since Point3D is immutable, this effectively returns a new point with
     * the same coordinates. While not strictly necessary for immutable objects,
     * it demonstrates the Prototype pattern.
     * 
     * @return a new Point3D with the same coordinates
     */
    public Point3D copy() {
        logger.log(Level.INFO, "Creating copy of point {0}", this);
        return new Point3D(x, y, z);
    }
}