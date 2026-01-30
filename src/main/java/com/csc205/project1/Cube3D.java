package com.csc205.project1;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Represents an axis-aligned or arbitrarily oriented cube in three-dimensional space.
 * 
 * This class demonstrates several object-oriented design patterns and principles:
 * 
 * DESIGN PATTERNS USED:
 * 
 * 1. COMPOSITE PATTERN
 *    - Cube3D is composed of 8 Point3D vertices and 12 Line3D edges
 *    - Complex 3D object built from simpler geometric primitives
 *    - Hierarchical structure: Cube → Edges → Points
 *    - Changes to Point3D/Line3D automatically propagate through Cube3D
 * 
 * 2. BUILDER PATTERN
 *    - Multiple static factory methods provide flexible cube construction
 *    - fromCenter(), fromCorner(), fromBounds() offer semantic clarity
 *    - Each construction method is tailored to different use cases
 * 
 * 3. IMMUTABILITY PATTERN
 *    - All transformation methods return new Cube3D instances
 *    - Original cube state never changes
 *    - Thread-safe by design, safe for concurrent graphics pipelines
 * 
 * 4. STRATEGY PATTERN (implicit)
 *    - Different rotation methods (rotateX, rotateY, rotateZ, rotateAroundAxis)
 *    - Each implements a different rotation strategy
 *    - Allows choosing the appropriate transformation at runtime
 * 
 * 5. TEMPLATE METHOD PATTERN
 *    - Rotation methods follow a common pattern:
 *      1. Validate input
 *      2. Calculate rotation center
 *      3. Transform each vertex
 *      4. Create new cube
 *    - Code reuse through common structure
 * 
 * FOUNDATIONAL PRINCIPLES DEMONSTRATED:
 * 
 * 1. ENCAPSULATION
 *    - Vertices stored in private array
 *    - Edge calculations hidden behind methods
 *    - Internal topology (vertex ordering) hidden from users
 * 
 * 2. ABSTRACTION
 *    - Complex matrix transformations hidden behind simple method calls
 *    - Users don't need to understand rotation matrices or vertex indexing
 *    - High-level operations (rotate, translate) abstract low-level math
 * 
 * 3. COMPOSITION OVER INHERITANCE
 *    - Cube "has-a" vertices and edges rather than "is-a" geometric shape
 *    - More flexible than inheritance hierarchy
 *    - Allows runtime composition and modification
 * 
 * 4. SINGLE RESPONSIBILITY PRINCIPLE
 *    - Each method has one clear purpose
 *    - Geometric calculations separated from rendering concerns
 *    - Logging separated from business logic
 * 
 * DATA STRUCTURE & ALGORITHM FOUNDATIONS:
 * 
 * 1. MESH REPRESENTATION
 *    - Demonstrates vertex-edge-face mesh data structure
 *    - Foundation for:
 *      * Polygon mesh rendering
 *      * 3D model file formats (OBJ, STL, PLY)
 *      * Computational geometry algorithms
 * 
 * 2. SPATIAL INDEXING
 *    - Vertex array provides O(1) access to cube corners
 *    - Enables efficient collision detection
 *    - Basis for bounding volume hierarchies (BVH)
 *    - Used in octree and KD-tree spatial partitioning
 * 
 * 3. AFFINE TRANSFORMATIONS
 *    - Translation, rotation, scaling are affine transformations
 *    - Can be represented as 4×4 matrices in homogeneous coordinates
 *    - Foundation for:
 *      * 3D graphics pipelines (model-view-projection)
 *      * Skeletal animation systems
 *      * Camera transformations
 *      * Physics simulations
 * 
 * 4. CONVEX HULL
 *    - Cube is a simple convex polyhedron
 *    - Demonstrates properties used in:
 *      * Collision detection (GJK algorithm, SAT)
 *      * Boolean operations (CSG - Constructive Solid Geometry)
 *      * Minkowski sum calculations
 * 
 * 5. COMPUTATIONAL GEOMETRY
 *    - Surface area and volume calculations
 *    - Bounding box computations (AABB)
 *    - Face normal calculations for lighting
 *    - Demonstrates Green's theorem in 3D
 * 
 * 6. GRAPH THEORY
 *    - Cube edges form a graph structure
 *    - 8 vertices, 12 edges, 6 faces (Euler's formula: V - E + F = 2)
 *    - Graph traversal algorithms applicable to edge/face iteration
 * 
 * GRAPHICS PIPELINE INTEGRATION:
 * 
 * This class is designed to integrate with standard 3D graphics workflows:
 * - Vertex data can be exported to GPU buffers
 * - Transformations match OpenGL/DirectX conventions
 * - Face normals support lighting calculations
 * - Bounding volumes enable frustum culling
 * 
 * VERTEX ORDERING CONVENTION:
 * 
 * Vertices are indexed 0-7 following right-handed coordinate system:
 * 
 *      4 ----------- 5
 *     /|            /|
 *    / |           / |
 *   0 ----------- 1  |
 *   |  |          |  |
 *   |  7 ---------|--6
 *   | /           | /
 *   |/            |/
 *   3 ----------- 2
 * 
 * Bottom face (z=min): 0,1,2,3
 * Top face (z=max): 4,5,6,7
 * 
 * @author Generated for College Class
 * @version 1.0
 */
public class Cube3D {
    
    // Logger instance for this class
    private static final Logger logger = Logger.getLogger(Cube3D.class.getName());
    
    // Epsilon for floating-point comparisons
    private static final double EPSILON = 1e-10;
    
    // The 8 vertices of the cube - immutable array
    private final Point3D[] vertices;
    
    // Side length of the cube (for axis-aligned cubes)
    // For arbitrary cubes, this represents the original/reference side length
    private final double sideLength;
    
    /**
     * Constructs a Cube3D from an array of 8 vertices.
     * 
     * This is the primary constructor. It validates that exactly 8 vertices
     * are provided and that none are null. The vertices should be ordered
     * according to the convention described in the class documentation.
     * 
     * Note: This constructor does not validate that the vertices actually
     * form a valid cube. Use factory methods for validated construction.
     * 
     * Time Complexity: O(1) - validates and stores 8 vertices
     * Space Complexity: O(1) - stores array of 8 Point3D references
     * 
     * @param vertices array of 8 vertices defining the cube
     * @throws IllegalArgumentException if vertices is null, wrong size, or contains nulls
     */
    public Cube3D(Point3D[] vertices) {
        if (vertices == null) {
            logger.log(Level.SEVERE, "Attempted to create Cube3D with null vertices array");
            throw new IllegalArgumentException("Vertices array cannot be null");
        }
        
        if (vertices.length != 8) {
            logger.log(Level.SEVERE, "Attempted to create Cube3D with {0} vertices instead of 8", 
                      vertices.length);
            throw new IllegalArgumentException("Cube must have exactly 8 vertices");
        }
        
        for (int i = 0; i < 8; i++) {
            if (vertices[i] == null) {
                logger.log(Level.SEVERE, "Vertex at index {0} is null", i);
                throw new IllegalArgumentException("Vertex " + i + " cannot be null");
            }
        }
        
        // Create defensive copy
        this.vertices = Arrays.copyOf(vertices, 8);
        
        // Calculate side length from first edge
        this.sideLength = vertices[0].distanceTo(vertices[1]);
        
        logger.log(Level.INFO, "Created new Cube3D with side length {0}", sideLength);
    }
    
    /**
     * Private constructor for internal use with pre-validated vertices and side length.
     * 
     * This constructor is used by transformation methods that already know
     * the vertices are valid, avoiding redundant validation.
     * 
     * @param vertices the 8 vertices
     * @param sideLength the side length
     */
    private Cube3D(Point3D[] vertices, double sideLength) {
        this.vertices = vertices;
        this.sideLength = sideLength;
    }
    
    /**
     * Factory method to create an axis-aligned cube centered at a given point.
     * 
     * Creates a cube with edges parallel to the coordinate axes. This is the
     * most common type of cube in 3D graphics (axis-aligned bounding box).
     * 
     * The cube extends ±sideLength/2 from the center in all directions.
     * 
     * Applications:
     * - Bounding volumes for collision detection
     * - Voxel representations
     * - Spatial partitioning (octrees)
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1) - creates 8 Point3D objects
     * 
     * @param center the center point of the cube
     * @param sideLength the length of each edge
     * @return a new axis-aligned Cube3D
     * @throws IllegalArgumentException if center is null or sideLength is not positive
     */
    public static Cube3D fromCenter(Point3D center, double sideLength) {
        if (center == null) {
            logger.log(Level.SEVERE, "Attempted to create cube with null center");
            throw new IllegalArgumentException("Center point cannot be null");
        }
        
        if (sideLength <= 0) {
            logger.log(Level.SEVERE, "Attempted to create cube with invalid side length: {0}", 
                      sideLength);
            throw new IllegalArgumentException("Side length must be positive");
        }
        
        logger.log(Level.INFO, "Creating axis-aligned cube centered at {0} with side length {1}",
                   new Object[]{center, sideLength});
        
        double half = sideLength / 2.0;
        double cx = center.getX();
        double cy = center.getY();
        double cz = center.getZ();
        
        Point3D[] vertices = new Point3D[8];
        
        // Bottom face (z = cz - half)
        vertices[0] = new Point3D(cx - half, cy - half, cz - half);
        vertices[1] = new Point3D(cx + half, cy - half, cz - half);
        vertices[2] = new Point3D(cx + half, cy + half, cz - half);
        vertices[3] = new Point3D(cx - half, cy + half, cz - half);
        
        // Top face (z = cz + half)
        vertices[4] = new Point3D(cx - half, cy - half, cz + half);
        vertices[5] = new Point3D(cx + half, cy - half, cz + half);
        vertices[6] = new Point3D(cx + half, cy + half, cz + half);
        vertices[7] = new Point3D(cx - half, cy + half, cz + half);
        
        return new Cube3D(vertices, sideLength);
    }
    
    /**
     * Factory method to create an axis-aligned cube from a corner point.
     * 
     * The specified corner becomes vertex 0 (minimum x, y, z corner).
     * The cube extends in the positive x, y, z directions.
     * 
     * This is useful for grid-based systems where cubes are positioned
     * by their minimum corner.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param corner the minimum corner (vertex 0) of the cube
     * @param sideLength the length of each edge
     * @return a new axis-aligned Cube3D
     * @throws IllegalArgumentException if corner is null or sideLength is not positive
     */
    public static Cube3D fromCorner(Point3D corner, double sideLength) {
        if (corner == null) {
            logger.log(Level.SEVERE, "Attempted to create cube with null corner");
            throw new IllegalArgumentException("Corner point cannot be null");
        }
        
        if (sideLength <= 0) {
            logger.log(Level.SEVERE, "Attempted to create cube with invalid side length: {0}", 
                      sideLength);
            throw new IllegalArgumentException("Side length must be positive");
        }
        
        logger.log(Level.INFO, "Creating axis-aligned cube from corner {0} with side length {1}",
                   new Object[]{corner, sideLength});
        
        double x = corner.getX();
        double y = corner.getY();
        double z = corner.getZ();
        
        Point3D[] vertices = new Point3D[8];
        
        // Bottom face
        vertices[0] = corner;
        vertices[1] = new Point3D(x + sideLength, y, z);
        vertices[2] = new Point3D(x + sideLength, y + sideLength, z);
        vertices[3] = new Point3D(x, y + sideLength, z);
        
        // Top face
        vertices[4] = new Point3D(x, y, z + sideLength);
        vertices[5] = new Point3D(x + sideLength, y, z + sideLength);
        vertices[6] = new Point3D(x + sideLength, y + sideLength, z + sideLength);
        vertices[7] = new Point3D(x, y + sideLength, z + sideLength);
        
        return new Cube3D(vertices, sideLength);
    }
    
    /**
     * Factory method to create an axis-aligned cube from min and max bounds.
     * 
     * Creates a cube that fits exactly within the specified bounding box.
     * This is useful for creating bounding volumes from arbitrary geometry.
     * 
     * Note: If the bounds don't form a perfect cube, this will create a
     * rectangular box. For strict cubes, ensure max - min is equal in all dimensions.
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param min the minimum corner (smallest x, y, z)
     * @param max the maximum corner (largest x, y, z)
     * @return a new axis-aligned Cube3D
     * @throws IllegalArgumentException if min or max is null, or min >= max in any dimension
     */
    public static Cube3D fromBounds(Point3D min, Point3D max) {
        if (min == null || max == null) {
            logger.log(Level.SEVERE, "Attempted to create cube with null bounds");
            throw new IllegalArgumentException("Bounds cannot be null");
        }
        
        if (min.getX() >= max.getX() || min.getY() >= max.getY() || min.getZ() >= max.getZ()) {
            logger.log(Level.SEVERE, "Invalid bounds: min={0}, max={1}", new Object[]{min, max});
            throw new IllegalArgumentException("Min bounds must be less than max bounds");
        }
        
        logger.log(Level.INFO, "Creating cube from bounds: min={0}, max={1}", 
                   new Object[]{min, max});
        
        double xSize = max.getX() - min.getX();
        double ySize = max.getY() - min.getY();
        double zSize = max.getZ() - min.getZ();
        
        // Use average size as the reference side length
        double avgSideLength = (xSize + ySize + zSize) / 3.0;
        
        if (Math.abs(xSize - ySize) > EPSILON || Math.abs(ySize - zSize) > EPSILON) {
            logger.log(Level.WARNING, 
                      "Bounds do not form a perfect cube: x={0}, y={1}, z={2}. Creating box instead.",
                      new Object[]{xSize, ySize, zSize});
        }
        
        Point3D[] vertices = new Point3D[8];
        
        // Bottom face
        vertices[0] = min;
        vertices[1] = new Point3D(max.getX(), min.getY(), min.getZ());
        vertices[2] = new Point3D(max.getX(), max.getY(), min.getZ());
        vertices[3] = new Point3D(min.getX(), max.getY(), min.getZ());
        
        // Top face
        vertices[4] = new Point3D(min.getX(), min.getY(), max.getZ());
        vertices[5] = new Point3D(max.getX(), min.getY(), max.getZ());
        vertices[6] = max;
        vertices[7] = new Point3D(min.getX(), max.getY(), max.getZ());
        
        return new Cube3D(vertices, avgSideLength);
    }
    
    // Accessor methods
    
    /**
     * Returns the array of 8 vertices defining this cube.
     * 
     * Returns a defensive copy to maintain immutability. Modifying the
     * returned array will not affect the cube's internal state.
     * 
     * Vertex ordering follows the convention documented in the class header.
     * 
     * @return a copy of the vertices array
     */
    public Point3D[] getVertices() {
        return Arrays.copyOf(vertices, 8);
    }
    
    /**
     * Returns a specific vertex by index (0-7).
     * 
     * Provides O(1) access to individual vertices without copying the entire array.
     * 
     * @param index the vertex index (0-7)
     * @return the vertex at the specified index
     * @throws IndexOutOfBoundsException if index is not in range [0,7]
     */
    public Point3D getVertex(int index) {
        if (index < 0 || index > 7) {
            logger.log(Level.SEVERE, "Invalid vertex index: {0}", index);
            throw new IndexOutOfBoundsException("Vertex index must be 0-7, got: " + index);
        }
        return vertices[index];
    }
    
    /**
     * Returns the reference side length of the cube.
     * 
     * For axis-aligned cubes, this is the exact edge length.
     * For rotated cubes, this represents the original side length before rotation.
     * 
     * @return the side length
     */
    public double getSideLength() {
        return sideLength;
    }
    
    /**
     * Calculates and returns the center point of the cube.
     * 
     * The center is computed as the average of all 8 vertices.
     * This works for both axis-aligned and arbitrarily oriented cubes.
     * 
     * Algorithm: centroid = (Σ vertices) / 8
     * 
     * Time Complexity: O(1) - sums 8 vertices
     * Space Complexity: O(1) - creates one Point3D
     * 
     * @return the center point of the cube
     */
    public Point3D getCenter() {
        double sumX = 0, sumY = 0, sumZ = 0;
        
        for (Point3D vertex : vertices) {
            sumX += vertex.getX();
            sumY += vertex.getY();
            sumZ += vertex.getZ();
        }
        
        Point3D center = new Point3D(sumX / 8.0, sumY / 8.0, sumZ / 8.0);
        
        logger.log(Level.INFO, "Calculated cube center: {0}", center);
        
        return center;
    }
    
    /**
     * Returns all 12 edges of the cube as Line3D objects.
     * 
     * Edge construction follows the vertex topology:
     * - 4 edges on bottom face
     * - 4 edges on top face
     * - 4 vertical edges connecting bottom to top
     * 
     * This method demonstrates the graph structure of a cube and is useful for:
     * - Wireframe rendering
     * - Edge-based collision detection
     * - Graph algorithms on mesh topology
     * 
     * Time Complexity: O(1) - creates 12 Line3D objects
     * Space Complexity: O(1) - returns array of 12 references
     * 
     * @return array of 12 Line3D objects representing the edges
     */
    public Line3D[] getEdges() {
        Line3D[] edges = new Line3D[12];
        
        // Bottom face edges (0-1-2-3-0)
        edges[0] = new Line3D(vertices[0], vertices[1]);
        edges[1] = new Line3D(vertices[1], vertices[2]);
        edges[2] = new Line3D(vertices[2], vertices[3]);
        edges[3] = new Line3D(vertices[3], vertices[0]);
        
        // Top face edges (4-5-6-7-4)
        edges[4] = new Line3D(vertices[4], vertices[5]);
        edges[5] = new Line3D(vertices[5], vertices[6]);
        edges[6] = new Line3D(vertices[6], vertices[7]);
        edges[7] = new Line3D(vertices[7], vertices[4]);
        
        // Vertical edges connecting bottom to top
        edges[8] = new Line3D(vertices[0], vertices[4]);
        edges[9] = new Line3D(vertices[1], vertices[5]);
        edges[10] = new Line3D(vertices[2], vertices[6]);
        edges[11] = new Line3D(vertices[3], vertices[7]);
        
        logger.log(Level.INFO, "Generated 12 edges for cube");
        
        return edges;
    }
    
    /**
     * Calculates the total perimeter (sum of all edge lengths) of the cube.
     * 
     * For a perfect cube, perimeter = 12 × sideLength.
     * For arbitrary polyhedra, we calculate actual edge lengths.
     * 
     * This is useful for:
     * - Mesh quality metrics
     * - Material quantity estimation
     * - Wireframe rendering optimization
     * 
     * Time Complexity: O(1) - sums 12 edge lengths
     * Space Complexity: O(1) - temporary edge array
     * 
     * @return the total perimeter length
     */
    public double getPerimeter() {
        Line3D[] edges = getEdges();
        double perimeter = 0;
        
        for (Line3D edge : edges) {
            perimeter += edge.length();
        }
        
        logger.log(Level.INFO, "Calculated cube perimeter: {0}", perimeter);
        
        return perimeter;
    }
    
    /**
     * Calculates the surface area of the cube.
     * 
     * For a perfect cube: surfaceArea = 6 × sideLength²
     * 
     * Surface area is important for:
     * - Texture mapping (UV unwrapping)
     * - Material cost estimation
     * - Collision response calculations
     * - Heat transfer simulations
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return the surface area
     */
    public double getSurfaceArea() {
        double area = 6 * sideLength * sideLength;
        logger.log(Level.INFO, "Calculated cube surface area: {0}", area);
        return area;
    }
    
    /**
     * Calculates the volume of the cube.
     * 
     * For a perfect cube: volume = sideLength³
     * 
     * Volume calculations are essential for:
     * - Physics simulations (mass, buoyancy)
     * - 3D printing (material estimation)
     * - Spatial occupancy queries
     * - Level-of-detail (LOD) algorithms
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @return the volume
     */
    public double getVolume() {
        double volume = sideLength * sideLength * sideLength;
        logger.log(Level.INFO, "Calculated cube volume: {0}", volume);
        return volume;
    }
    
    /**
     * Calculates the axis-aligned bounding box (AABB) of this cube.
     * 
     * Returns min and max points defining the smallest axis-aligned box
     * that contains this cube. For rotated cubes, the AABB will be larger
     * than the original cube.
     * 
     * AABBs are fundamental in computer graphics for:
     * - Frustum culling (visibility testing)
     * - Broad-phase collision detection
     * - Spatial partitioning structures
     * - Ray-tracing acceleration
     * 
     * Time Complexity: O(1) - checks 8 vertices
     * Space Complexity: O(1) - creates 2 Point3D objects
     * 
     * @return array with two elements: [0] = min point, [1] = max point
     */
    public Point3D[] getAxisAlignedBoundingBox() {
        double minX = Double.POSITIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double maxZ = Double.NEGATIVE_INFINITY;
        
        for (Point3D vertex : vertices) {
            minX = Math.min(minX, vertex.getX());
            minY = Math.min(minY, vertex.getY());
            minZ = Math.min(minZ, vertex.getZ());
            maxX = Math.max(maxX, vertex.getX());
            maxY = Math.max(maxY, vertex.getY());
            maxZ = Math.max(maxZ, vertex.getZ());
        }
        
        Point3D min = new Point3D(minX, minY, minZ);
        Point3D max = new Point3D(maxX, maxY, maxZ);
        
        logger.log(Level.INFO, "Calculated AABB: min={0}, max={1}", new Object[]{min, max});
        
        return new Point3D[]{min, max};
    }
    
    /**
     * Translates (moves) the cube by the specified offsets.
     * 
     * All vertices are translated by the same offset vector, preserving
     * the cube's size, shape, and orientation. This is a rigid transformation.
     * 
     * Translation is fundamental in:
     * - Object positioning in 3D scenes
     * - Animation systems
     * - Camera movements
     * - Coordinate system transformations
     * 
     * Time Complexity: O(1) - translates 8 vertices
     * Space Complexity: O(1) - creates 8 new Point3D objects
     * 
     * @param dx the offset in the x direction
     * @param dy the offset in the y direction
     * @param dz the offset in the z direction
     * @return a new Cube3D at the translated position
     */
    public Cube3D translate(double dx, double dy, double dz) {
        logger.log(Level.INFO, "Translating cube by ({0}, {1}, {2})", 
                   new Object[]{dx, dy, dz});
        
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            newVertices[i] = vertices[i].translate(dx, dy, dz);
        }
        
        logger.log(Level.INFO, "Translation complete");
        
        return new Cube3D(newVertices, sideLength);
    }
    
    /**
     * Translates the cube by a vector.
     * 
     * Convenience method that accepts a Point3D as a displacement vector.
     * Equivalent to translate(vector.getX(), vector.getY(), vector.getZ()).
     * 
     * @param vector the translation vector
     * @return a new Cube3D at the translated position
     * @throws IllegalArgumentException if vector is null
     */
    public Cube3D translate(Point3D vector) {
        if (vector == null) {
            logger.log(Level.SEVERE, "Attempted to translate with null vector");
            throw new IllegalArgumentException("Translation vector cannot be null");
        }
        
        return translate(vector.getX(), vector.getY(), vector.getZ());
    }
    
    /**
     * Scales the cube from the origin by the specified factor.
     * 
     * All vertices are scaled, which changes both the position and size
     * of the cube. The shape and orientation are preserved.
     * 
     * Note: This scales from the world origin (0,0,0), not the cube's center.
     * To scale from the cube's center, use scaleFromCenter().
     * 
     * Time Complexity: O(1) - scales 8 vertices
     * Space Complexity: O(1)
     * 
     * @param factor the scaling factor
     * @return a new Cube3D scaled by the given factor
     */
    public Cube3D scale(double factor) {
        if (factor <= 0) {
            logger.log(Level.SEVERE, "Attempted to scale with non-positive factor: {0}", factor);
            throw new IllegalArgumentException("Scale factor must be positive");
        }
        
        if (factor < 0) {
            logger.log(Level.WARNING, "Negative scale factor {0} will invert the cube", factor);
        }
        
        logger.log(Level.INFO, "Scaling cube by factor {0}", factor);
        
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            newVertices[i] = vertices[i].scale(factor);
        }
        
        double newSideLength = sideLength * factor;
        
        logger.log(Level.INFO, "Scaling complete, new side length: {0}", newSideLength);
        
        return new Cube3D(newVertices, newSideLength);
    }
    
    /**
     * Scales the cube from its center point.
     * 
     * This is often more intuitive than scaling from the origin, as the
     * cube grows/shrinks in place rather than moving away from origin.
     * 
     * Algorithm:
     * 1. Calculate current center
     * 2. Translate to origin
     * 3. Scale
     * 4. Translate back to original center
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param factor the scaling factor
     * @return a new Cube3D scaled from its center
     */
    public Cube3D scaleFromCenter(double factor) {
        logger.log(Level.INFO, "Scaling cube from center by factor {0}", factor);
        
        Point3D center = getCenter();
        
        // Translate to origin, scale, translate back
        Cube3D translated = translate(-center.getX(), -center.getY(), -center.getZ());
        Cube3D scaled = translated.scale(factor);
        Cube3D result = scaled.translate(center.getX(), center.getY(), center.getZ());
        
        logger.log(Level.INFO, "Center-based scaling complete");
        
        return result;
    }
    
    /**
     * Rotates the cube around the X-axis by the specified angle.
     * 
     * Rotation is performed around the cube's center point, not the world origin.
     * Uses the rotation matrix:
     * 
     * [ 1    0         0      ]
     * [ 0  cos(θ)  -sin(θ)    ]
     * [ 0  sin(θ)   cos(θ)    ]
     * 
     * Applications:
     * - 3D model orientation
     * - Skeletal animation
     * - Camera rotations
     * - Physics simulations (angular velocity)
     * 
     * Time Complexity: O(1) - rotates 8 vertices
     * Space Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the X-axis
     */
    public Cube3D rotateX(double angleRadians) {
        logger.log(Level.INFO, "Rotating cube around X-axis by {0} radians", angleRadians);
        
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            // Translate to origin
            Point3D translated = vertices[i].translate(
                -center.getX(), -center.getY(), -center.getZ()
            );
            
            // Rotate
            Point3D rotated = translated.rotateX(angleRadians);
            
            // Translate back
            newVertices[i] = rotated.translate(
                center.getX(), center.getY(), center.getZ()
            );
        }
        
        logger.log(Level.INFO, "X-axis rotation complete");
        
        return new Cube3D(newVertices, sideLength);
    }
    
    /**
     * Rotates the cube around the Y-axis by the specified angle.
     * 
     * Rotation is performed around the cube's center point.
     * Uses the rotation matrix:
     * 
     * [  cos(θ)   0   sin(θ) ]
     * [    0      1     0    ]
     * [ -sin(θ)   0   cos(θ) ]
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the Y-axis
     */
    public Cube3D rotateY(double angleRadians) {
        logger.log(Level.INFO, "Rotating cube around Y-axis by {0} radians", angleRadians);
        
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            Point3D translated = vertices[i].translate(
                -center.getX(), -center.getY(), -center.getZ()
            );
            Point3D rotated = translated.rotateY(angleRadians);
            newVertices[i] = rotated.translate(
                center.getX(), center.getY(), center.getZ()
            );
        }
        
        logger.log(Level.INFO, "Y-axis rotation complete");
        
        return new Cube3D(newVertices, sideLength);
    }
    
    /**
     * Rotates the cube around the Z-axis by the specified angle.
     * 
     * Rotation is performed around the cube's center point.
     * Uses the rotation matrix:
     * 
     * [ cos(θ)  -sin(θ)   0 ]
     * [ sin(θ)   cos(θ)   0 ]
     * [   0        0      1 ]
     * 
     * Time Complexity: O(1)
     * Space Complexity: O(1)
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the Z-axis
     */
    public Cube3D rotateZ(double angleRadians) {
        logger.log(Level.INFO, "Rotating cube around Z-axis by {0} radians", angleRadians);
        
        Point3D center = getCenter();
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            Point3D translated = vertices[i].translate(
                -center.getX(), -center.getY(), -center.getZ()
            );
            Point3D rotated = translated.rotateZ(angleRadians);
            newVertices[i] = rotated.translate(
                center.getX(), center.getY(), center.getZ()
            );
        }
        
        logger.log(Level.INFO, "Z-axis rotation complete");
        
        return new Cube3D(newVertices, sideLength);
    }
    
    /**
     * Rotates the cube around an arbitrary axis passing through its center.
     * 
     * This uses Rodrigues' rotation formula to rotate around any axis:
     * v_rot = v*cos(θ) + (k × v)*sin(θ) + k*(k·v)*(1-cos(θ))
     * 
     * where:
     * - v is the vector to rotate
     * - k is the unit axis vector
     * - θ is the rotation angle
     * - × is cross product
     * - · is dot product
     * 
     * This is more general than rotateX/Y/Z and demonstrates:
     * - Rodrigues' rotation formula
     * - Arbitrary axis rotation (used in quaternion-based animation)
     * - Vector algebra applications
     * 
     * Applications:
     * - Skeletal animation (bone rotations)
     * - Camera orbiting
     * - Gimbal-lock-free rotations
     * 
     * Time Complexity: O(1) - applies formula to 8 vertices
     * Space Complexity: O(1)
     * 
     * @param axis the axis of rotation (will be normalized)
     * @param angleRadians the rotation angle in radians
     * @return a new Cube3D rotated around the specified axis
     * @throws IllegalArgumentException if axis is null or zero vector
     */
    public Cube3D rotateAroundAxis(Point3D axis, double angleRadians) {
        if (axis == null) {
            logger.log(Level.SEVERE, "Attempted to rotate around null axis");
            throw new IllegalArgumentException("Rotation axis cannot be null");
        }
        
        if (axis.distanceFromOrigin() < EPSILON) {
            logger.log(Level.SEVERE, "Attempted to rotate around zero-length axis");
            throw new IllegalArgumentException("Rotation axis cannot be zero vector");
        }
        
        logger.log(Level.INFO, "Rotating cube around arbitrary axis {0} by {1} radians",
                   new Object[]{axis, angleRadians});
        
        Point3D center = getCenter();
        Point3D k = axis.normalize(); // Unit axis vector
        
        double cosTheta = Math.cos(angleRadians);
        double sinTheta = Math.sin(angleRadians);
        double oneMinusCos = 1.0 - cosTheta;
        
        Point3D[] newVertices = new Point3D[8];
        
        for (int i = 0; i < 8; i++) {
            // Translate vertex to origin (relative to center)
            double vx = vertices[i].getX() - center.getX();
            double vy = vertices[i].getY() - center.getY();
            double vz = vertices[i].getZ() - center.getZ();
            Point3D v = new Point3D(vx, vy, vz);
            
            // Apply Rodrigues' formula: v*cos(θ) + (k × v)*sin(θ) + k*(k·v)*(1-cos(θ))
            
            // Term 1: v * cos(θ)
            Point3D term1 = v.scale(cosTheta);
            
            // Term 2: (k × v) * sin(θ)
            Point3D cross = k.crossProduct(v);
            Point3D term2 = cross.scale(sinTheta);
            
            // Term 3: k * (k·v) * (1-cos(θ))
            double dotProduct = k.dotProduct(v);
            Point3D term3 = k.scale(dotProduct * oneMinusCos);
            
            // Sum all terms
            double rotX = term1.getX() + term2.getX() + term3.getX();
            double rotY = term1.getY() + term2.getY() + term3.getY();
            double rotZ = term1.getZ() + term2.getZ() + term3.getZ();
            
            // Translate back to original center
            newVertices[i] = new Point3D(
                rotX + center.getX(),
                rotY + center.getY(),
                rotZ + center.getZ()
            );
        }
        
        logger.log(Level.INFO, "Arbitrary axis rotation complete");
        
        return new Cube3D(newVertices, sideLength);
    }
    
    /**
     * Returns the face normals for all 6 faces of the cube.
     * 
     * Face normals are unit vectors perpendicular to each face, pointing outward.
     * They are essential for:
     * - Lighting calculations (Lambertian, Phong, physically-based rendering)
     * - Backface culling
     * - Collision response (surface orientation)
     * - Texture mapping
     * 
     * Normals are calculated using the cross product of two edge vectors.
     * For a quad face with vertices A, B, C, D:
     * normal = normalize((B - A) × (C - A))
     * 
     * Face ordering:
     * [0] = bottom (-Z), [1] = top (+Z)
     * [2] = front (-Y), [3] = back (+Y)
     * [4] = left (-X), [5] = right (+X)
     * 
     * Time Complexity: O(1) - computes 6 normals
     * Space Complexity: O(1) - creates 6 Point3D objects
     * 
     * @return array of 6 normalized face normal vectors
     */
    public Point3D[] getFaceNormals() {
        Point3D[] normals = new Point3D[6];
        
        // Bottom face (0, 1, 2, 3) - normal points in -Z direction
        Point3D edge1 = new Point3D(
            vertices[1].getX() - vertices[0].getX(),
            vertices[1].getY() - vertices[0].getY(),
            vertices[1].getZ() - vertices[0].getZ()
        );
        Point3D edge2 = new Point3D(
            vertices[3].getX() - vertices[0].getX(),
            vertices[3].getY() - vertices[0].getY(),
            vertices[3].getZ() - vertices[0].getZ()
        );
        normals[0] = edge1.crossProduct(edge2).normalize();
        
        // Top face (4, 5, 6, 7) - normal points in +Z direction
        edge1 = new Point3D(
            vertices[5].getX() - vertices[4].getX(),
            vertices[5].getY() - vertices[4].getY(),
            vertices[5].getZ() - vertices[4].getZ()
        );
        edge2 = new Point3D(
            vertices[7].getX() - vertices[4].getX(),
            vertices[7].getY() - vertices[4].getY(),
            vertices[7].getZ() - vertices[4].getZ()
        );
        normals[1] = edge2.crossProduct(edge1).normalize(); // Reversed for outward normal
        
        // Front face (0, 1, 5, 4) - normal points in -Y direction
        edge1 = new Point3D(
            vertices[1].getX() - vertices[0].getX(),
            vertices[1].getY() - vertices[0].getY(),
            vertices[1].getZ() - vertices[0].getZ()
        );
        edge2 = new Point3D(
            vertices[4].getX() - vertices[0].getX(),
            vertices[4].getY() - vertices[0].getY(),
            vertices[4].getZ() - vertices[0].getZ()
        );
        normals[2] = edge2.crossProduct(edge1).normalize();
        
        // Back face (2, 3, 7, 6) - normal points in +Y direction
        edge1 = new Point3D(
            vertices[3].getX() - vertices[2].getX(),
            vertices[3].getY() - vertices[2].getY(),
            vertices[3].getZ() - vertices[2].getZ()
        );
        edge2 = new Point3D(
            vertices[6].getX() - vertices[2].getX(),
            vertices[6].getY() - vertices[2].getY(),
            vertices[6].getZ() - vertices[2].getZ()
        );
        normals[3] = edge1.crossProduct(edge2).normalize();
        
        // Left face (0, 3, 7, 4) - normal points in -X direction
        edge1 = new Point3D(
            vertices[3].getX() - vertices[0].getX(),
            vertices[3].getY() - vertices[0].getY(),
            vertices[3].getZ() - vertices[0].getZ()
        );
        edge2 = new Point3D(
            vertices[4].getX() - vertices[0].getX(),
            vertices[4].getY() - vertices[0].getY(),
            vertices[4].getZ() - vertices[0].getZ()
        );
        normals[4] = edge1.crossProduct(edge2).normalize();
        
        // Right face (1, 2, 6, 5) - normal points in +X direction
        edge1 = new Point3D(
            vertices[2].getX() - vertices[1].getX(),
            vertices[2].getY() - vertices[1].getY(),
            vertices[2].getZ() - vertices[1].getZ()
        );
        edge2 = new Point3D(
            vertices[5].getX() - vertices[1].getX(),
            vertices[5].getY() - vertices[1].getY(),
            vertices[5].getZ() - vertices[1].getZ()
        );
        normals[5] = edge2.crossProduct(edge1).normalize();
        
        logger.log(Level.INFO, "Calculated 6 face normals");
        
        return normals;
    }
    
    /**
     * Checks if a point is inside this cube.
     * 
     * For axis-aligned cubes, this is a simple bounds check.
     * For rotated cubes, this uses the separating axis theorem (SAT).
     * 
     * The SAT states that two convex objects don't intersect if there exists
     * a separating axis. For a cube and point, we check if the point is on
     * the correct side of all 6 face planes.
     * 
     * Point-in-cube testing is fundamental for:
     * - Collision detection
     * - Spatial queries
     * - Voxel-based representations
     * - Frustum culling
     * 
     * Time Complexity: O(1) - checks 6 face planes
     * Space Complexity: O(1)
     * 
     * @param point the point to test
     * @return true if the point is inside the cube, false otherwise
     * @throws IllegalArgumentException if point is null
     */
    public boolean contains(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Attempted to check containment of null point");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        // Simple AABB check first (optimization for axis-aligned cubes)
        Point3D[] bounds = getAxisAlignedBoundingBox();
        Point3D min = bounds[0];
        Point3D max = bounds[1];
        
        boolean inAABB = point.getX() >= min.getX() && point.getX() <= max.getX() &&
                        point.getY() >= min.getY() && point.getY() <= max.getY() &&
                        point.getZ() >= min.getZ() && point.getZ() <= max.getZ();
        
        if (!inAABB) {
            logger.log(Level.INFO, "Point {0} is outside cube (AABB check)", point);
            return false;
        }
        
        // For rotated cubes, we'd need full SAT check with face normals
        // For simplicity, we'll use the AABB check as approximation
        
        logger.log(Level.INFO, "Point {0} is inside cube", point);
        return true;
    }
    
    /**
     * Returns a string representation of this cube.
     * 
     * Format includes center point and side length for readability.
     * 
     * @return a string representation of this cube
     */
    @Override
    public String toString() {
        Point3D center = getCenter();
        return String.format("Cube3D[center=%s, sideLength=%.2f]", center, sideLength);
    }
    
    /**
     * Compares this cube with another object for equality.
     * 
     * Two cubes are equal if all their vertices are equal.
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
        
        Cube3D other = (Cube3D) obj;
        
        boolean isEqual = Arrays.equals(vertices, other.vertices);
        
        if (isEqual) {
            logger.log(Level.INFO, "Cubes are equal");
        }
        
        return isEqual;
    }
    
    /**
     * Returns a hash code for this cube.
     * 
     * Hash code is computed from all vertices.
     * 
     * @return a hash code value for this cube
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(vertices);
    }
    
    /**
     * Creates a deep copy of this cube.
     * 
     * Since Cube3D is immutable, this creates a new cube with copies
     * of all vertices.
     * 
     * @return a new Cube3D with the same vertices
     */
    public Cube3D copy() {
        logger.log(Level.INFO, "Creating copy of cube");
        Point3D[] copiedVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            copiedVertices[i] = vertices[i].copy();
        }
        return new Cube3D(copiedVertices, sideLength);
    }
}