package com.csc205.project1;

import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a line segment in three-dimensional Euclidean space.
 * 
 * This class demonstrates several object-oriented design patterns and principles:
 * 
 * DESIGN PATTERNS USED:
 * 
 * 1. COMPOSITION PATTERN
 *    - Line3D is composed of two Point3D objects (start and end points)
 *    - Demonstrates "has-a" relationship rather than "is-a"
 *    - Promotes code reuse and maintains single responsibility
 *    - Changes to Point3D behavior automatically propagate to Line3D
 * 
 * 2. VALUE OBJECT PATTERN
 *    - Immutable state (final fields) ensures thread safety
 *    - Lines are compared by value, not identity
 *    - Safe to use in collections and as hash map keys
 * 
 * 3. FACADE PATTERN
 *    - Complex geometric calculations are hidden behind simple method interfaces
 *    - Users don't need to understand vector mathematics, parametric equations,
 *      or cross products to use the class
 * 
 * 4. FACTORY PATTERN (static factory methods)
 *    - Static creation methods provide semantic clarity
 *    - Alternative construction paths (fromPoints, fromDirection, etc.)
 *    - Can return cached instances or perform validation before construction
 * 
 * FOUNDATIONAL PRINCIPLES DEMONSTRATED:
 * 
 * 1. ENCAPSULATION
 *    - Private fields with public accessor methods
 *    - Internal mathematical complexity hidden from users
 *    - Logging concern separated via Logger
 * 
 * 2. ABSTRACTION
 *    - Complex algorithms (skew line distance, parametric intersection)
 *      abstracted into simple method calls
 *    - Mathematical details hidden behind intuitive interfaces
 * 
 * 3. DELEGATION
 *    - Line3D delegates point-specific operations to Point3D class
 *    - Demonstrates proper use of composition
 *    - Avoids code duplication
 * 
 * 4. IMMUTABILITY
 *    - All operations return new Line3D instances
 *    - Thread-safe by design
 *    - Simplifies reasoning about state
 * 
 * DATA STRUCTURE & ALGORITHM FOUNDATIONS:
 * 
 * 1. GEOMETRIC PRIMITIVES
 *    - Lines are fundamental for computational geometry algorithms:
 *      * Line sweep algorithms
 *      * Convex hull construction
 *      * Visibility graphs
 *      * Ray tracing and rendering
 * 
 * 2. SPATIAL RELATIONSHIPS
 *    - Line-line intersection testing is basis for:
 *      * Collision detection in games/physics
 *      * Boolean operations on polygons
 *      * Network topology analysis
 * 
 * 3. PARAMETRIC REPRESENTATION
 *    - Lines represented as P(t) = P0 + t * direction
 *    - Foundation for:
 *      * Bezier curves and splines
 *      * Animation paths
 *      * Interpolation algorithms
 * 
 * 4. DISTANCE METRICS
 *    - Shortest distance between skew lines uses:
 *      * Vector projections
 *      * Cross products
 *      * Linear algebra fundamentals
 *    - Critical for proximity queries in spatial databases
 * 
 * 5. COMPUTATIONAL GEOMETRY
 *    - Methods demonstrate classic algorithms:
 *      * Line-line distance (closest point of approach)
 *      * Point-line distance (perpendicular projection)
 *      * Parallel/perpendicular testing
 * 
 * @author Generated for College Class
 * @version 1.0
 */
public class Line3D {
    
    // Logger instance for this class - demonstrates Singleton pattern
    private static final Logger logger = Logger.getLogger(Line3D.class.getName());
    
    // Epsilon for floating-point comparisons
    private static final double EPSILON = 1e-10;
    
    // Endpoints of the line segment - immutable for thread safety
    private final Point3D start;
    private final Point3D end;
    
    /**
     * Constructs a new Line3D with the specified start and end points.
     * 
     * This is the primary constructor. It validates that the two points are
     * distinct (not the same point), as a line requires two different points.
     * 
     * The composition pattern is demonstrated here: Line3D "has-a" relationship
     * with Point3D, rather than inheriting from it.
     * 
     * Time Complexity: O(1) - constant time initialization
     * Space Complexity: O(1) - stores references to two Point3D objects
     * 
     * @param start the starting point of the line segment
     * @param end the ending point of the line segment
     * @throws IllegalArgumentException if either point is null or points are identical
     */
    public Line3D(Point3D start, Point3D end) {
        if (start == null || end == null) {
            logger.log(Level.SEVERE, "Attempted to create Line3D with null point(s)");
            throw new IllegalArgumentException("Line endpoints cannot be null");
        }
        
        if (start.equals(end)) {
            logger.log(Level.SEVERE, "Attempted to create Line3D with identical points: {0}", start);
            throw new IllegalArgumentException("Line endpoints must be distinct points");
        }
        
        this.start = start;
        this.end = end;
        
        logger.log(Level.INFO, "Created new Line3D from {0} to {1}", 
                   new Object[]{start, end});
    }
    
    /**
     * Factory method to create a line from two points.
     * 
     * This static factory method provides a more semantic way to construct
     * a line, making the code more readable: Line3D.fromPoints(p1, p2)
     * 
     * Factory methods offer advantages over constructors:
     * - Descriptive names clarify intent
     * - Can perform additional validation or preprocessing
     * - Can return cached instances (though we don't here)
     * 
     * @param p1 the first point
     * @param p2 the second point
     * @return a new Line3D from p1 to p2
     */
    public static Line3D fromPoints(Point3D p1, Point3D p2) {
        logger.log(Level.INFO, "Creating line via factory method fromPoints");
        return new Line3D(p1, p2);
    }
    
    /**
     * Factory method to create a line from a starting point and direction vector.
     * 
     * Creates a line segment of specified length starting at the origin point
     * and extending in the given direction. The direction vector doesn't need
     * to be normalized.
     * 
     * This demonstrates the parametric line representation:
     * P(t) = origin + t * direction, where t = length
     * 
     * @param origin the starting point
     * @param direction the direction vector (as a Point3D)
     * @param length the length of the line segment
     * @return a new Line3D from origin extending in direction
     * @throws IllegalArgumentException if inputs are invalid
     */
    public static Line3D fromDirection(Point3D origin, Point3D direction, double length) {
        if (origin == null || direction == null) {
            logger.log(Level.SEVERE, "Null point provided to fromDirection factory");
            throw new IllegalArgumentException("Origin and direction cannot be null");
        }
        
        if (length <= 0) {
            logger.log(Level.SEVERE, "Invalid length {0} provided to fromDirection", length);
            throw new IllegalArgumentException("Length must be positive");
        }
        
        logger.log(Level.INFO, "Creating line from direction: origin={0}, direction={1}, length={2}",
                   new Object[]{origin, direction, length});
        
        // Normalize direction and scale by length
        Point3D normalizedDir = direction.normalize();
        Point3D scaledDir = normalizedDir.scale(length);
        Point3D endpoint = origin.translate(scaledDir.getX(), scaledDir.getY(), scaledDir.getZ());
        
        return new Line3D(origin, endpoint);
    }
    
    // Accessor methods - demonstrate encapsulation
    
    /**
     * Returns the starting point of this line segment.
     * 
     * Returns a reference to the immutable Point3D, which is safe because
     * Point3D itself is immutable. If Point3D were mutable, we would need
     * to return a defensive copy.
     * 
     * @return the start point
     */
    public Point3D getStart() {
        return start;
    }
    
    /**
     * Returns the ending point of this line segment.
     * 
     * @return the end point
     */
    public Point3D getEnd() {
        return end;
    }
    
    /**
     * Calculates and returns the length of this line segment.
     * 
     * Length is computed using the Euclidean distance between endpoints.
     * This delegates to Point3D's distanceTo() method, demonstrating
     * proper use of composition and avoiding code duplication.
     * 
     * Time Complexity: O(1) - constant number of operations
     * Space Complexity: O(1) - only uses primitive return value
     * 
     * @return the length of the line segment
     */
    public double length() {
        double len = start.distanceTo(end);
        logger.log(Level.INFO, "Calculated line length: {0}", len);
        return len;
    }
    
    /**
     * Returns the direction vector of this line.
     * 
     * The direction vector points from start to end and represents the
     * line's orientation in 3D space. This is not normalized, so its
     * magnitude equals the line's length.
     * 
     * Direction vectors are fundamental in:
     * - Ray tracing algorithms
     * - Physics simulations (velocity, force)
     * - Computer graphics transformations
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1) - creates one Point3D
     * 
     * @return a Point3D representing the direction vector
     */
    public Point3D getDirectionVector() {
        double dx = end.getX() - start.getX();
        double dy = end.getY() - start.getY();
        double dz = end.getZ() - start.getZ();
        
        Point3D direction = new Point3D(dx, dy, dz);
        logger.log(Level.INFO, "Calculated direction vector: {0}", direction);
        
        return direction;
    }
    
    /**
     * Returns the normalized (unit length) direction vector of this line.
     * 
     * A unit direction vector has magnitude 1 and indicates only the
     * direction, not the magnitude. This is useful when you only care
     * about orientation, not length.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return a Point3D representing the normalized direction vector
     */
    public Point3D getUnitDirectionVector() {
        Point3D direction = getDirectionVector();
        Point3D normalized = direction.normalize();
        logger.log(Level.INFO, "Calculated unit direction vector: {0}", normalized);
        return normalized;
    }
    
    /**
     * Calculates the midpoint of this line segment.
     * 
     * The midpoint is equidistant from both endpoints and lies on the line.
     * Midpoint calculation is used in:
     * - Binary space partitioning
     * - Subdivision surfaces
     * - Collision detection (bounding volume hierarchies)
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return a Point3D at the midpoint of the line
     */
    public Point3D midpoint() {
        Point3D mid = start.midpoint(end);
        logger.log(Level.INFO, "Calculated line midpoint: {0}", mid);
        return mid;
    }
    
    /**
     * Calculates a point at a given parameter t along the line.
     * 
     * Uses the parametric line equation: P(t) = start + t * (end - start)
     * - t = 0 returns the start point
     * - t = 1 returns the end point
     * - 0 < t < 1 returns a point between start and end
     * - t < 0 or t > 1 returns a point on the infinite line extension
     * 
     * Parametric equations are fundamental in:
     * - Animation and interpolation
     * - Bezier curves and splines
     * - Path planning algorithms
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param t the parameter value
     * @return the point at parameter t
     */
    public Point3D pointAt(double t) {
        if (t < 0 || t > 1) {
            logger.log(Level.WARNING, "Parameter t={0} is outside [0,1], point will be off segment", t);
        }
        
        double x = start.getX() + t * (end.getX() - start.getX());
        double y = start.getY() + t * (end.getY() - start.getY());
        double z = start.getZ() + t * (end.getZ() - start.getZ());
        
        Point3D point = new Point3D(x, y, z);
        logger.log(Level.INFO, "Calculated point at t={0}: {1}", new Object[]{t, point});
        
        return point;
    }
    
    /**
     * Calculates the shortest distance from a point to this line segment.
     * 
     * This uses vector projection to find the perpendicular distance.
     * 
     * Algorithm:
     * 1. Project the point onto the infinite line
     * 2. Clamp the projection to the segment [0, 1]
     * 3. Calculate distance from point to closest point on segment
     * 
     * Applications:
     * - Collision detection (point-line proximity)
     * - Nearest neighbor queries
     * - Pathfinding (distance to obstacles)
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param point the point to measure distance from
     * @return the shortest distance from the point to this line segment
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to calculate distance to null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        Point3D direction = getDirectionVector();
        double lengthSquared = Math.pow(start.distanceTo(end), 2);
        
        // Vector from start to point
        double dx = point.getX() - start.getX();
        double dy = point.getY() - start.getY();
        double dz = point.getZ() - start.getZ();
        Point3D startToPoint = new Point3D(dx, dy, dz);
        
        // Project onto line: t = (point - start) · direction / |direction|²
        double t = startToPoint.dotProduct(direction) / lengthSquared;
        
        // Clamp t to [0, 1] to stay on segment
        t = Math.max(0, Math.min(1, t));
        
        // Find closest point on segment
        Point3D closestPoint = pointAt(t);
        
        double distance = point.distanceTo(closestPoint);
        
        logger.log(Level.INFO, "Distance from point {0} to line: {1}", 
                   new Object[]{point, distance});
        
        return distance;
    }
    
    /**
     * Calculates the shortest distance between this line and another line.
     * 
     * This is one of the most complex geometric calculations, handling three cases:
     * 1. PARALLEL LINES: Distance between any point on one line to the other
     * 2. INTERSECTING LINES: Distance is zero
     * 3. SKEW LINES: Uses the formula: |((P2-P1) · (d1 × d2))| / |d1 × d2|
     * 
     * The algorithm demonstrates:
     * - Vector cross products
     * - Dot products
     * - Handling of degenerate cases (parallel lines)
     * - Numerical stability considerations
     * 
     * Applications:
     * - Collision detection between moving objects
     * - Robot path planning (clearance checking)
     * - Computer graphics (closest approach)
     * - Network topology (minimum spanning trees)
     * 
     * Time Complexity: O(1) - fixed number of vector operations
     * Space Complexity: O(1) - creates temporary Point3D objects
     * 
     * Mathematical Background:
     * For two lines L1 = P1 + s*d1 and L2 = P2 + t*d2:
     * - If d1 × d2 = 0, lines are parallel
     * - Otherwise, distance = |(P2-P1) · (d1 × d2)| / |d1 × d2|
     * 
     * @param other the other line to calculate distance to
     * @return the shortest distance between the two lines
     * @throws IllegalArgumentException if other is null
     */
    public double shortestDistanceTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to calculate distance to null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        logger.log(Level.INFO, "Calculating shortest distance between lines");
        
        // Get direction vectors
        Point3D d1 = this.getDirectionVector();
        Point3D d2 = other.getDirectionVector();
        
        // Vector between start points
        double dx = other.start.getX() - this.start.getX();
        double dy = other.start.getY() - this.start.getY();
        double dz = other.start.getZ() - this.start.getZ();
        Point3D w = new Point3D(dx, dy, dz);
        
        // Cross product of direction vectors
        Point3D crossProduct = d1.crossProduct(d2);
        double crossMagnitude = crossProduct.distanceFromOrigin();
        
        // Check if lines are parallel (cross product ≈ 0)
        if (crossMagnitude < EPSILON) {
            logger.log(Level.INFO, "Lines are parallel, using point-to-line distance");
            
            // Lines are parallel - distance from any point on one line to the other
            double dist = this.distanceToPoint(other.start);
            
            logger.log(Level.INFO, "Distance between parallel lines: {0}", dist);
            return dist;
        }
        
        // Lines are skew - use the formula: |(w · (d1 × d2))| / |d1 × d2|
        double numerator = Math.abs(w.dotProduct(crossProduct));
        double distance = numerator / crossMagnitude;
        
        logger.log(Level.INFO, "Distance between skew lines: {0}", distance);
        
        return distance;
    }
    
    /**
     * Determines if this line is parallel to another line.
     * 
     * Two lines are parallel if their direction vectors are parallel,
     * which occurs when their cross product is zero (or very close to zero).
     * 
     * Uses epsilon comparison for floating-point robustness.
     * 
     * Parallel line detection is important for:
     * - Geometric constraint solving
     * - CAD/CAM applications
     * - Collision detection optimization
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param other the other line to check
     * @return true if the lines are parallel, false otherwise
     * @throws IllegalArgumentException if other is null
     */
    public boolean isParallelTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to check parallelism with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        Point3D d1 = this.getDirectionVector();
        Point3D d2 = other.getDirectionVector();
        
        Point3D cross = d1.crossProduct(d2);
        double crossMagnitude = cross.distanceFromOrigin();
        
        boolean parallel = crossMagnitude < EPSILON;
        
        logger.log(Level.INFO, "Lines are parallel: {0}", parallel);
        
        return parallel;
    }
    
    /**
     * Determines if this line is perpendicular to another line.
     * 
     * Two lines are perpendicular if their direction vectors are perpendicular,
     * which occurs when their dot product is zero (or very close to zero).
     * 
     * Uses epsilon comparison for floating-point robustness.
     * 
     * Perpendicularity testing is used in:
     * - Orthogonal coordinate system verification
     * - Right angle detection in geometry
     * - Surface normal calculations
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param other the other line to check
     * @return true if the lines are perpendicular, false otherwise
     * @throws IllegalArgumentException if other is null
     */
    public boolean isPerpendicularTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Attempted to check perpendicularity with null line");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        Point3D d1 = this.getDirectionVector();
        Point3D d2 = other.getDirectionVector();
        
        double dotProduct = d1.dotProduct(d2);
        boolean perpendicular = Math.abs(dotProduct) < EPSILON;
        
        logger.log(Level.INFO, "Lines are perpendicular: {0}", perpendicular);
        
        return perpendicular;
    }
    
    /**
     * Translates (moves) this line by the specified offsets.
     * 
     * Both endpoints are translated by the same offset, preserving the
     * line's length and direction. This is a rigid transformation.
     * 
     * Translation is fundamental in:
     * - Animation and motion
     * - Coordinate system transformations
     * - Object positioning in 3D scenes
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param dx the offset in the x direction
     * @param dy the offset in the y direction
     * @param dz the offset in the z direction
     * @return a new Line3D at the translated position
     */
    public Line3D translate(double dx, double dy, double dz) {
        logger.log(Level.INFO, "Translating line by ({0}, {1}, {2})", 
                   new Object[]{dx, dy, dz});
        
        Point3D newStart = start.translate(dx, dy, dz);
        Point3D newEnd = end.translate(dx, dy, dz);
        
        Line3D translated = new Line3D(newStart, newEnd);
        
        logger.log(Level.INFO, "Translation complete");
        
        return translated;
    }
    
    /**
     * Scales this line from the origin by the specified factor.
     * 
     * Both endpoints are scaled, which changes both the position and length
     * of the line. The direction is preserved (unless factor is negative).
     * 
     * Scaling transformations are used in:
     * - Zoom operations
     * - Model scaling in 3D graphics
     * - Coordinate system conversions
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param factor the scaling factor
     * @return a new Line3D scaled by the given factor
     */
    public Line3D scale(double factor) {
        if (factor < 0) {
            logger.log(Level.WARNING, "Scaling with negative factor {0} will invert the line", 
                      factor);
        }
        
        logger.log(Level.INFO, "Scaling line by factor {0}", factor);
        
        Point3D newStart = start.scale(factor);
        Point3D newEnd = end.scale(factor);
        
        Line3D scaled = new Line3D(newStart, newEnd);
        
        logger.log(Level.INFO, "Scaling complete, new length: {0}", scaled.length());
        
        return scaled;
    }
    
    /**
     * Reverses the direction of this line.
     * 
     * Creates a new line with start and end points swapped. The line segment
     * is geometrically the same but has opposite direction.
     * 
     * Direction reversal is used in:
     * - Graph algorithms (reversing edges)
     * - Path reversal in navigation
     * - Normal flipping in graphics
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return a new Line3D with reversed direction
     */
    public Line3D reverse() {
        logger.log(Level.INFO, "Reversing line direction");
        Line3D reversed = new Line3D(end, start);
        logger.log(Level.INFO, "Line reversed: {0} → {1}", new Object[]{start, end});
        return reversed;
    }
    
    /**
     * Extends this line segment by the specified factor.
     * 
     * Factor > 1: Makes the line longer
     * Factor = 1: No change
     * 0 < Factor < 1: Makes the line shorter
     * 
     * The line is extended symmetrically from both ends, maintaining
     * the same center point and direction.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param factor the extension factor
     * @return a new Line3D extended by the given factor
     * @throws IllegalArgumentException if factor is not positive
     */
    public Line3D extend(double factor) {
        if (factor <= 0) {
            logger.log(Level.SEVERE, "Invalid extension factor: {0}", factor);
            throw new IllegalArgumentException("Extension factor must be positive");
        }
        
        logger.log(Level.INFO, "Extending line by factor {0}", factor);
        
        Point3D mid = midpoint();
        Point3D direction = getDirectionVector();
        
        // Scale direction by half the extension factor (symmetric extension)
        double halfExtension = factor / 2.0;
        Point3D scaledDir = direction.scale(halfExtension);
        
        Point3D newStart = new Point3D(
            mid.getX() - scaledDir.getX(),
            mid.getY() - scaledDir.getY(),
            mid.getZ() - scaledDir.getZ()
        );
        
        Point3D newEnd = new Point3D(
            mid.getX() + scaledDir.getX(),
            mid.getY() + scaledDir.getY(),
            mid.getZ() + scaledDir.getZ()
        );
        
        Line3D extended = new Line3D(newStart, newEnd);
        
        logger.log(Level.INFO, "Extension complete, new length: {0}", extended.length());
        
        return extended;
    }
    
    /**
     * Returns a string representation of this line.
     * 
     * Format: "Line3D[start → end]"
     * 
     * String representation is important for:
     * - Debugging and logging
     * - User interface display
     * - Serialization for persistence
     * 
     * @return a string representation of this line
     */
    @Override
    public String toString() {
        return String.format("Line3D[%s → %s]", start, end);
    }
    
    /**
     * Compares this line with another object for equality.
     * 
     * Two lines are equal if they have the same start and end points.
     * Note: This checks directed equality. Line(A,B) ≠ Line(B,A).
     * 
     * For undirected equality (where direction doesn't matter), you would
     * need to check both (start==other.start && end==other.end) OR
     * (start==other.end && end==other.start).
     * 
     * Proper equals() implementation is critical for:
     * - Using objects in HashMaps and HashSets
     * - Graph edge comparison in algorithms
     * - Duplicate detection
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
        
        Line3D other = (Line3D) obj;
        
        boolean isEqual = start.equals(other.start) && end.equals(other.end);
        
        if (isEqual) {
            logger.log(Level.INFO, "Lines are equal: {0} == {1}", 
                      new Object[]{this, other});
        }
        
        return isEqual;
    }
    
    /**
     * Returns a hash code for this line.
     * 
     * Hash codes are essential for using objects in hash-based collections.
     * The contract states that equal objects must have equal hash codes.
     * 
     * This implementation combines the hash codes of both endpoints.
     * 
     * @return a hash code value for this line
     */
    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + start.hashCode();
        result = 31 * result + end.hashCode();
        return result;
    }
    
    /**
     * Creates a deep copy of this line.
     * 
     * Since Line3D is immutable (and composed of immutable Point3D objects),
     * this effectively returns a new line with the same endpoints.
     * 
     * Demonstrates the Prototype pattern for object cloning.
     * 
     * @return a new Line3D with the same endpoints
     */
    public Line3D copy() {
        logger.log(Level.INFO, "Creating copy of line {0}", this);
        return new Line3D(start.copy(), end.copy());
    }
}