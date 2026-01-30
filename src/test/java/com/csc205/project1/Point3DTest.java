package com.csc205.project1;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive unit test suite for the Point3D class.
 *
 * This test class demonstrates best practices in unit testing:
 * - Organized into nested test classes for logical grouping
 * - Tests both normal cases and edge cases
 * - Uses descriptive test names and display names
 * - Follows AAA pattern (Arrange, Act, Assert)
 * - Tests boundary conditions and special values
 * - Verifies exception handling
 * - Tests immutability guarantees
 *
 * Test Coverage:
 * - Constructors and factory methods
 * - Accessor methods (getters)
 * - Distance calculations (Euclidean, Manhattan, from origin)
 * - Rotation operations (around each axis)
 * - Translation and scaling
 * - Vector operations (dot product, cross product, normalize)
 * - Utility methods (midpoint, equals, hashCode, toString)
 * - Edge cases (zero vectors, null inputs, extreme values)
 *
 * @author Generated for College Class
 * @version 1.0
 */
class Point3DTest {

        // Tolerance for floating-point comparisons
        private static final double EPSILON = 1e-9;

        // Common test points used across multiple tests
        private Point3D origin;
        private Point3D unitX;
        private Point3D unitY;
        private Point3D unitZ;
        private Point3D point123;
        private Point3D pointNegative;

        /**
         * Sets up common test fixtures before each test.
         *
         * This demonstrates the DRY (Don't Repeat Yourself) principle
         * and ensures consistent test data across all tests.
         */
        @BeforeEach
        void setUp() {
            origin = new Point3D(0, 0, 0);
            unitX = new Point3D(1, 0, 0);
            unitY = new Point3D(0, 1, 0);
            unitZ = new Point3D(0, 0, 1);
            point123 = new Point3D(1, 2, 3);
            pointNegative = new Point3D(-1, -2, -3);
        }

        /**
         * Tests for constructor and factory methods.
         */
        @Nested
        @DisplayName("Constructor and Factory Method Tests")
        class ConstructorTests {

            @Test
            @DisplayName("Should create point with specified coordinates")
            void testConstructorWithCoordinates() {
                Point3D point = new Point3D(1.5, 2.5, 3.5);

                assertEquals(1.5, point.getX(), EPSILON);
                assertEquals(2.5, point.getY(), EPSILON);
                assertEquals(3.5, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should create point at origin with default constructor")
            void testDefaultConstructor() {
                Point3D point = new Point3D();

                assertEquals(0.0, point.getX(), EPSILON);
                assertEquals(0.0, point.getY(), EPSILON);
                assertEquals(0.0, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle negative coordinates")
            void testConstructorWithNegativeCoordinates() {
                Point3D point = new Point3D(-5.5, -10.2, -15.8);

                assertEquals(-5.5, point.getX(), EPSILON);
                assertEquals(-10.2, point.getY(), EPSILON);
                assertEquals(-15.8, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle zero coordinates")
            void testConstructorWithZeros() {
                Point3D point = new Point3D(0, 0, 0);

                assertEquals(0.0, point.getX(), EPSILON);
                assertEquals(0.0, point.getY(), EPSILON);
                assertEquals(0.0, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should create origin point via factory method")
            void testOriginFactoryMethod() {
                Point3D point = Point3D.origin();

                assertEquals(0.0, point.getX(), EPSILON);
                assertEquals(0.0, point.getY(), EPSILON);
                assertEquals(0.0, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle very large coordinates")
            void testConstructorWithLargeValues() {
                Point3D point = new Point3D(1e10, 1e10, 1e10);

                assertEquals(1e10, point.getX(), EPSILON);
                assertEquals(1e10, point.getY(), EPSILON);
                assertEquals(1e10, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle very small coordinates")
            void testConstructorWithSmallValues() {
                Point3D point = new Point3D(1e-10, 1e-10, 1e-10);

                assertEquals(1e-10, point.getX(), EPSILON);
                assertEquals(1e-10, point.getY(), EPSILON);
                assertEquals(1e-10, point.getZ(), EPSILON);
            }
        }

        /**
         * Tests for distance calculation methods.
         */
        @Nested
        @DisplayName("Distance Calculation Tests")
        class DistanceTests {

            @Test
            @DisplayName("Should calculate distance between two points correctly")
            void testDistanceTo() {
                Point3D p1 = new Point3D(0, 0, 0);
                Point3D p2 = new Point3D(3, 4, 0);

                double distance = p1.distanceTo(p2);

                assertEquals(5.0, distance, EPSILON);
            }

            @Test
            @DisplayName("Should calculate 3D distance correctly")
            void testDistanceTo3D() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(4, 6, 8);

                // Distance = sqrt((4-1)² + (6-2)² + (8-3)²) = sqrt(9 + 16 + 25) = sqrt(50)
                double expected = Math.sqrt(50);
                double distance = p1.distanceTo(p2);

                assertEquals(expected, distance, EPSILON);
            }

            @Test
            @DisplayName("Should return zero distance for same point")
            void testDistanceToSamePoint() {
                double distance = point123.distanceTo(point123);

                assertEquals(0.0, distance, EPSILON);
            }

            @Test
            @DisplayName("Should be symmetric (d(A,B) = d(B,A))")
            void testDistanceSymmetry() {
                double d1 = point123.distanceTo(pointNegative);
                double d2 = pointNegative.distanceTo(point123);

                assertEquals(d1, d2, EPSILON);
            }

            @Test
            @DisplayName("Should throw exception for null point in distanceTo")
            void testDistanceToNull() {
                assertThrows(IllegalArgumentException.class, () -> {
                    point123.distanceTo(null);
                });
            }

            @Test
            @DisplayName("Should calculate Manhattan distance correctly")
            void testManhattanDistance() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(4, 6, 8);

                // Manhattan = |4-1| + |6-2| + |8-3| = 3 + 4 + 5 = 12
                double distance = p1.manhattanDistanceTo(p2);

                assertEquals(12.0, distance, EPSILON);
            }

            @Test
            @DisplayName("Should handle negative coordinates in Manhattan distance")
            void testManhattanDistanceWithNegatives() {
                Point3D p1 = new Point3D(-1, -2, -3);
                Point3D p2 = new Point3D(1, 2, 3);

                // Manhattan = |1-(-1)| + |2-(-2)| + |3-(-3)| = 2 + 4 + 6 = 12
                double distance = p1.manhattanDistanceTo(p2);

                assertEquals(12.0, distance, EPSILON);
            }

            @Test
            @DisplayName("Should throw exception for null point in manhattanDistanceTo")
            void testManhattanDistanceToNull() {
                assertThrows(IllegalArgumentException.class, () -> {
                    point123.manhattanDistanceTo(null);
                });
            }

            @Test
            @DisplayName("Should calculate distance from origin correctly")
            void testDistanceFromOrigin() {
                Point3D point = new Point3D(3, 4, 0);

                double distance = point.distanceFromOrigin();

                assertEquals(5.0, distance, EPSILON);
            }

            @Test
            @DisplayName("Should return zero for origin point distance from origin")
            void testOriginDistanceFromOrigin() {
                double distance = origin.distanceFromOrigin();

                assertEquals(0.0, distance, EPSILON);
            }

            @Test
            @DisplayName("Should calculate 3D distance from origin")
            void testDistanceFromOrigin3D() {
                Point3D point = new Point3D(1, 2, 2);

                // Distance = sqrt(1² + 2² + 2²) = sqrt(9) = 3
                double distance = point.distanceFromOrigin();

                assertEquals(3.0, distance, EPSILON);
            }
        }

        /**
         * Tests for rotation operations.
         */
        @Nested
        @DisplayName("Rotation Tests")
        class RotationTests {

            @Test
            @DisplayName("Should rotate point around X-axis by 90 degrees")
            void testRotateX90Degrees() {
                Point3D point = new Point3D(0, 1, 0);
                Point3D rotated = point.rotateX(Math.PI / 2);

                assertEquals(0.0, rotated.getX(), EPSILON);
                assertEquals(0.0, rotated.getY(), EPSILON);
                assertEquals(1.0, rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should rotate point around X-axis by 180 degrees")
            void testRotateX180Degrees() {
                Point3D point = new Point3D(0, 1, 1);
                Point3D rotated = point.rotateX(Math.PI);

                assertEquals(0.0, rotated.getX(), EPSILON);
                assertEquals(-1.0, rotated.getY(), EPSILON);
                assertEquals(-1.0, rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should not change X coordinate when rotating around X-axis")
            void testRotateXPreservesX() {
                Point3D rotated = point123.rotateX(Math.PI / 4);

                assertEquals(point123.getX(), rotated.getX(), EPSILON);
            }

            @Test
            @DisplayName("Should rotate point around Y-axis by 90 degrees")
            void testRotateY90Degrees() {
                Point3D point = new Point3D(1, 0, 0);
                Point3D rotated = point.rotateY(Math.PI / 2);

                assertEquals(0.0, rotated.getX(), EPSILON);
                assertEquals(0.0, rotated.getY(), EPSILON);
                assertEquals(-1.0, rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should not change Y coordinate when rotating around Y-axis")
            void testRotateYPreservesY() {
                Point3D rotated = point123.rotateY(Math.PI / 3);

                assertEquals(point123.getY(), rotated.getY(), EPSILON);
            }

            @Test
            @DisplayName("Should rotate point around Z-axis by 90 degrees")
            void testRotateZ90Degrees() {
                Point3D point = new Point3D(1, 0, 0);
                Point3D rotated = point.rotateZ(Math.PI / 2);

                assertEquals(0.0, rotated.getX(), EPSILON);
                assertEquals(1.0, rotated.getY(), EPSILON);
                assertEquals(0.0, rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should not change Z coordinate when rotating around Z-axis")
            void testRotateZPreservesZ() {
                Point3D rotated = point123.rotateZ(Math.PI / 6);

                assertEquals(point123.getZ(), rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should preserve distance from origin during rotation")
            void testRotationPreservesDistance() {
                double originalDistance = point123.distanceFromOrigin();

                Point3D rotatedX = point123.rotateX(Math.PI / 4);
                Point3D rotatedY = point123.rotateY(Math.PI / 4);
                Point3D rotatedZ = point123.rotateZ(Math.PI / 4);

                assertEquals(originalDistance, rotatedX.distanceFromOrigin(), EPSILON);
                assertEquals(originalDistance, rotatedY.distanceFromOrigin(), EPSILON);
                assertEquals(originalDistance, rotatedZ.distanceFromOrigin(), EPSILON);
            }

            @Test
            @DisplayName("Should handle rotation by zero radians")
            void testRotationByZero() {
                Point3D rotatedX = point123.rotateX(0);
                Point3D rotatedY = point123.rotateY(0);
                Point3D rotatedZ = point123.rotateZ(0);

                assertEquals(point123.getX(), rotatedX.getX(), EPSILON);
                assertEquals(point123.getY(), rotatedX.getY(), EPSILON);
                assertEquals(point123.getZ(), rotatedX.getZ(), EPSILON);

                assertEquals(point123.getX(), rotatedY.getX(), EPSILON);
                assertEquals(point123.getY(), rotatedY.getY(), EPSILON);
                assertEquals(point123.getZ(), rotatedY.getZ(), EPSILON);

                assertEquals(point123.getX(), rotatedZ.getX(), EPSILON);
                assertEquals(point123.getY(), rotatedZ.getY(), EPSILON);
                assertEquals(point123.getZ(), rotatedZ.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle rotation by 2π (full circle)")
            void testRotationByTwoPi() {
                Point3D rotated = point123.rotateX(2 * Math.PI);

                assertEquals(point123.getX(), rotated.getX(), EPSILON);
                assertEquals(point123.getY(), rotated.getY(), EPSILON);
                assertEquals(point123.getZ(), rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle negative rotation angles")
            void testNegativeRotation() {
                Point3D rotatedPositive = unitX.rotateZ(Math.PI / 2);
                Point3D rotatedNegative = unitX.rotateZ(-Math.PI / 2);

                // Rotating +90° and -90° should give opposite Y values
                assertEquals(rotatedPositive.getY(), -rotatedNegative.getY(), EPSILON);
            }
        }

        /**
         * Tests for translation operations.
         */
        @Nested
        @DisplayName("Translation Tests")
        class TranslationTests {

            @Test
            @DisplayName("Should translate point correctly")
            void testTranslate() {
                Point3D translated = point123.translate(1, 2, 3);

                assertEquals(2.0, translated.getX(), EPSILON);
                assertEquals(4.0, translated.getY(), EPSILON);
                assertEquals(6.0, translated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should translate by negative values")
            void testTranslateNegative() {
                Point3D translated = point123.translate(-1, -2, -3);

                assertEquals(0.0, translated.getX(), EPSILON);
                assertEquals(0.0, translated.getY(), EPSILON);
                assertEquals(0.0, translated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle translation by zero")
            void testTranslateByZero() {
                Point3D translated = point123.translate(0, 0, 0);

                assertEquals(point123.getX(), translated.getX(), EPSILON);
                assertEquals(point123.getY(), translated.getY(), EPSILON);
                assertEquals(point123.getZ(), translated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle large translation values")
            void testTranslateLargeValues() {
                Point3D translated = origin.translate(1e6, 1e6, 1e6);

                assertEquals(1e6, translated.getX(), EPSILON);
                assertEquals(1e6, translated.getY(), EPSILON);
                assertEquals(1e6, translated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Translation should be reversible")
            void testTranslationReversible() {
                Point3D translated = point123.translate(5, 10, 15);
                Point3D back = translated.translate(-5, -10, -15);

                assertEquals(point123.getX(), back.getX(), EPSILON);
                assertEquals(point123.getY(), back.getY(), EPSILON);
                assertEquals(point123.getZ(), back.getZ(), EPSILON);
            }
        }

        /**
         * Tests for scaling operations.
         */
        @Nested
        @DisplayName("Scaling Tests")
        class ScalingTests {

            @Test
            @DisplayName("Should scale point correctly")
            void testScale() {
                Point3D scaled = point123.scale(2);

                assertEquals(2.0, scaled.getX(), EPSILON);
                assertEquals(4.0, scaled.getY(), EPSILON);
                assertEquals(6.0, scaled.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should scale by fraction")
            void testScaleByFraction() {
                Point3D scaled = point123.scale(0.5);

                assertEquals(0.5, scaled.getX(), EPSILON);
                assertEquals(1.0, scaled.getY(), EPSILON);
                assertEquals(1.5, scaled.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle scale factor of 1")
            void testScaleByOne() {
                Point3D scaled = point123.scale(1);

                assertEquals(point123.getX(), scaled.getX(), EPSILON);
                assertEquals(point123.getY(), scaled.getY(), EPSILON);
                assertEquals(point123.getZ(), scaled.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle negative scale factor")
            void testScaleByNegative() {
                Point3D scaled = point123.scale(-1);

                assertEquals(-1.0, scaled.getX(), EPSILON);
                assertEquals(-2.0, scaled.getY(), EPSILON);
                assertEquals(-3.0, scaled.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should scale origin to origin")
            void testScaleOrigin() {
                Point3D scaled = origin.scale(5);

                assertEquals(0.0, scaled.getX(), EPSILON);
                assertEquals(0.0, scaled.getY(), EPSILON);
                assertEquals(0.0, scaled.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Scaling should preserve ratios")
            void testScalePreservesRatios() {
                Point3D scaled = point123.scale(3);

                // Ratios should be preserved: x:y:z = 1:2:3
                double ratio1 = scaled.getX() / scaled.getY();
                double ratio2 = point123.getX() / point123.getY();

                assertEquals(ratio2, ratio1, EPSILON);
            }
        }

        /**
         * Tests for vector operations.
         */
        @Nested
        @DisplayName("Vector Operation Tests")
        class VectorOperationTests {

            @Test
            @DisplayName("Should calculate dot product correctly")
            void testDotProduct() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(4, 5, 6);

                // Dot product = 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
                double dotProduct = p1.dotProduct(p2);

                assertEquals(32.0, dotProduct, EPSILON);
            }

            @Test
            @DisplayName("Dot product should be commutative")
            void testDotProductCommutative() {
                double d1 = point123.dotProduct(pointNegative);
                double d2 = pointNegative.dotProduct(point123);

                assertEquals(d1, d2, EPSILON);
            }

            @Test
            @DisplayName("Dot product of perpendicular vectors should be zero")
            void testDotProductPerpendicular() {
                double dotProduct = unitX.dotProduct(unitY);

                assertEquals(0.0, dotProduct, EPSILON);
            }

            @Test
            @DisplayName("Dot product with itself equals magnitude squared")
            void testDotProductWithSelf() {
                double dotProduct = point123.dotProduct(point123);
                double magnitudeSquared = Math.pow(point123.distanceFromOrigin(), 2);

                assertEquals(magnitudeSquared, dotProduct, EPSILON);
            }

            @Test
            @DisplayName("Should throw exception for null in dot product")
            void testDotProductNull() {
                assertThrows(IllegalArgumentException.class, () -> {
                    point123.dotProduct(null);
                });
            }

            @Test
            @DisplayName("Should calculate cross product correctly")
            void testCrossProduct() {
                Point3D result = unitX.crossProduct(unitY);

                assertEquals(0.0, result.getX(), EPSILON);
                assertEquals(0.0, result.getY(), EPSILON);
                assertEquals(1.0, result.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Cross product should be anti-commutative")
            void testCrossProductAntiCommutative() {
                Point3D cross1 = unitX.crossProduct(unitY);
                Point3D cross2 = unitY.crossProduct(unitX);

                assertEquals(cross1.getX(), -cross2.getX(), EPSILON);
                assertEquals(cross1.getY(), -cross2.getY(), EPSILON);
                assertEquals(cross1.getZ(), -cross2.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Cross product of parallel vectors should be zero")
            void testCrossProductParallel() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(2, 4, 6); // Parallel to p1

                Point3D cross = p1.crossProduct(p2);

                assertEquals(0.0, cross.getX(), EPSILON);
                assertEquals(0.0, cross.getY(), EPSILON);
                assertEquals(0.0, cross.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Cross product result should be perpendicular to both inputs")
            void testCrossProductPerpendicular() {
                Point3D cross = point123.crossProduct(pointNegative);

                // Cross product should be perpendicular (dot product = 0)
                double dot1 = cross.dotProduct(point123);
                double dot2 = cross.dotProduct(pointNegative);

                assertEquals(0.0, dot1, EPSILON);
                assertEquals(0.0, dot2, EPSILON);
            }

            @Test
            @DisplayName("Should throw exception for null in cross product")
            void testCrossProductNull() {
                assertThrows(IllegalArgumentException.class, () -> {
                    point123.crossProduct(null);
                });
            }

            @Test
            @DisplayName("Should normalize vector correctly")
            void testNormalize() {
                Point3D normalized = point123.normalize();

                double magnitude = normalized.distanceFromOrigin();

                assertEquals(1.0, magnitude, EPSILON);
            }

            @Test
            @DisplayName("Normalized vector should preserve direction")
            void testNormalizePreservesDirection() {
                Point3D normalized = point123.normalize();

                // Check that ratios are preserved (direction is same)
                double ratioOriginal = point123.getX() / point123.getY();
                double ratioNormalized = normalized.getX() / normalized.getY();

                assertEquals(ratioOriginal, ratioNormalized, EPSILON);
            }

            @Test
            @DisplayName("Should throw exception when normalizing zero vector")
            void testNormalizeZeroVector() {
                assertThrows(ArithmeticException.class, () -> {
                    origin.normalize();
                });
            }

            @Test
            @DisplayName("Normalizing a unit vector should return unit vector")
            void testNormalizeUnitVector() {
                Point3D normalized = unitX.normalize();

                assertEquals(1.0, normalized.getX(), EPSILON);
                assertEquals(0.0, normalized.getY(), EPSILON);
                assertEquals(0.0, normalized.getZ(), EPSILON);
            }
        }

        /**
         * Tests for utility methods.
         */
        @Nested
        @DisplayName("Utility Method Tests")
        class UtilityTests {

            @Test
            @DisplayName("Should calculate midpoint correctly")
            void testMidpoint() {
                Point3D p1 = new Point3D(0, 0, 0);
                Point3D p2 = new Point3D(4, 6, 8);

                Point3D mid = p1.midpoint(p2);

                assertEquals(2.0, mid.getX(), EPSILON);
                assertEquals(3.0, mid.getY(), EPSILON);
                assertEquals(4.0, mid.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Midpoint should be commutative")
            void testMidpointCommutative() {
                Point3D mid1 = point123.midpoint(pointNegative);
                Point3D mid2 = pointNegative.midpoint(point123);

                assertEquals(mid1.getX(), mid2.getX(), EPSILON);
                assertEquals(mid1.getY(), mid2.getY(), EPSILON);
                assertEquals(mid1.getZ(), mid2.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Midpoint with negative coordinates")
            void testMidpointWithNegatives() {
                Point3D p1 = new Point3D(-2, -4, -6);
                Point3D p2 = new Point3D(2, 4, 6);

                Point3D mid = p1.midpoint(p2);

                assertEquals(0.0, mid.getX(), EPSILON);
                assertEquals(0.0, mid.getY(), EPSILON);
                assertEquals(0.0, mid.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should throw exception for null in midpoint")
            void testMidpointNull() {
                assertThrows(IllegalArgumentException.class, () -> {
                    point123.midpoint(null);
                });
            }

            @Test
            @DisplayName("Should create correct copy")
            void testCopy() {
                Point3D copy = point123.copy();

                assertEquals(point123.getX(), copy.getX(), EPSILON);
                assertEquals(point123.getY(), copy.getY(), EPSILON);
                assertEquals(point123.getZ(), copy.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Copy should be equal but not same object")
            void testCopyNotSameObject() {
                Point3D copy = point123.copy();

                assertTrue(point123.equals(copy));
                assertNotSame(point123, copy);
            }
        }

        /**
         * Tests for equals and hashCode methods.
         */
        @Nested
        @DisplayName("Equals and HashCode Tests")
        class EqualsHashCodeTests {

            @Test
            @DisplayName("Should be equal to itself")
            void testEqualsReflexive() {
                assertTrue(point123.equals(point123));
            }

            @Test
            @DisplayName("Should be equal to point with same coordinates")
            void testEqualsWithSameCoordinates() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(1, 2, 3);

                assertTrue(p1.equals(p2));
            }

            @Test
            @DisplayName("Should be symmetric")
            void testEqualsSymmetric() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(1, 2, 3);

                assertTrue(p1.equals(p2));
                assertTrue(p2.equals(p1));
            }

            @Test
            @DisplayName("Should be transitive")
            void testEqualsTransitive() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(1, 2, 3);
                Point3D p3 = new Point3D(1, 2, 3);

                assertTrue(p1.equals(p2));
                assertTrue(p2.equals(p3));
                assertTrue(p1.equals(p3));
            }

            @Test
            @DisplayName("Should not be equal to null")
            void testEqualsNull() {
                assertFalse(point123.equals(null));
            }

            @Test
            @DisplayName("Should not be equal to different class")
            void testEqualsDifferentClass() {
                assertFalse(point123.equals("not a point"));
            }

            @Test
            @DisplayName("Should not be equal to point with different coordinates")
            void testNotEqualsDifferentCoordinates() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(1, 2, 4);

                assertFalse(p1.equals(p2));
            }

            @Test
            @DisplayName("Equal points should have same hash code")
            void testHashCodeConsistency() {
                Point3D p1 = new Point3D(1, 2, 3);
                Point3D p2 = new Point3D(1, 2, 3);

                assertEquals(p1.hashCode(), p2.hashCode());
            }

            @Test
            @DisplayName("Hash code should be consistent across calls")
            void testHashCodeStability() {
                int hash1 = point123.hashCode();
                int hash2 = point123.hashCode();

                assertEquals(hash1, hash2);
            }

            @Test
            @DisplayName("Different points should likely have different hash codes")
            void testHashCodeDifferent() {
                int hash1 = point123.hashCode();
                int hash2 = pointNegative.hashCode();

                assertNotEquals(hash1, hash2);
            }
        }

        /**
         * Tests for toString method.
         */
        @Nested
        @DisplayName("String Representation Tests")
        class ToStringTests {

            @Test
            @DisplayName("toString should contain all coordinates")
            void testToString() {
                String str = point123.toString();

                assertTrue(str.contains("1"));
                assertTrue(str.contains("2"));
                assertTrue(str.contains("3"));
            }

            @Test
            @DisplayName("toString should be non-null and non-empty")
            void testToStringNotEmpty() {
                String str = point123.toString();

                assertNotNull(str);
                assertFalse(str.isEmpty());
            }

            @Test
            @DisplayName("toString should handle negative coordinates")
            void testToStringNegative() {
                String str = pointNegative.toString();

                assertTrue(str.contains("-1") || str.contains("−1"));
                assertTrue(str.contains("-2") || str.contains("−2"));
                assertTrue(str.contains("-3") || str.contains("−3"));
            }
        }

        /**
         * Tests for immutability guarantees.
         */
        @Nested
        @DisplayName("Immutability Tests")
        class ImmutabilityTests {

            @Test
            @DisplayName("Rotation should not modify original point")
            void testRotationImmutability() {
                Point3D original = new Point3D(1, 2, 3);
                original.rotateX(Math.PI / 2);

                assertEquals(1.0, original.getX(), EPSILON);
                assertEquals(2.0, original.getY(), EPSILON);
                assertEquals(3.0, original.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Translation should not modify original point")
            void testTranslationImmutability() {
                Point3D original = new Point3D(1, 2, 3);
                original.translate(5, 5, 5);

                assertEquals(1.0, original.getX(), EPSILON);
                assertEquals(2.0, original.getY(), EPSILON);
                assertEquals(3.0, original.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Scaling should not modify original point")
            void testScalingImmutability() {
                Point3D original = new Point3D(1, 2, 3);
                original.scale(10);

                assertEquals(1.0, original.getX(), EPSILON);
                assertEquals(2.0, original.getY(), EPSILON);
                assertEquals(3.0, original.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Normalize should not modify original point")
            void testNormalizeImmutability() {
                Point3D original = new Point3D(3, 4, 0);
                double originalMagnitude = original.distanceFromOrigin();

                original.normalize();

                assertEquals(originalMagnitude, original.distanceFromOrigin(), EPSILON);
            }
        }

        /**
         * Tests for edge cases and boundary conditions.
         */
        @Nested
        @DisplayName("Edge Case Tests")
        class EdgeCaseTests {

            @Test
            @DisplayName("Should handle very small floating point values")
            void testVerySmallValues() {
                Point3D point = new Point3D(1e-15, 1e-15, 1e-15);

                assertNotNull(point);
                assertTrue(point.distanceFromOrigin() > 0);
            }

            @Test
            @DisplayName("Should handle very large floating point values")
            void testVeryLargeValues() {
                Point3D point = new Point3D(1e15, 1e15, 1e15);

                assertNotNull(point);
                assertTrue(point.distanceFromOrigin() > 0);
            }

            @Test
            @DisplayName("Should handle mixed positive and negative coordinates")
            void testMixedSigns() {
                Point3D point = new Point3D(-1, 2, -3);

                assertEquals(-1.0, point.getX(), EPSILON);
                assertEquals(2.0, point.getY(), EPSILON);
                assertEquals(-3.0, point.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle rotation of origin")
            void testRotateOrigin() {
                Point3D rotated = origin.rotateX(Math.PI / 2);

                assertEquals(0.0, rotated.getX(), EPSILON);
                assertEquals(0.0, rotated.getY(), EPSILON);
                assertEquals(0.0, rotated.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle cross product of identical vectors")
            void testCrossProductIdentical() {
                Point3D cross = point123.crossProduct(point123);

                assertEquals(0.0, cross.getX(), EPSILON);
                assertEquals(0.0, cross.getY(), EPSILON);
                assertEquals(0.0, cross.getZ(), EPSILON);
            }

            @Test
            @DisplayName("Should handle dot product with zero vector")
            void testDotProductWithZero() {
                double dot = point123.dotProduct(origin);

                assertEquals(0.0, dot, EPSILON);
            }
        }
    }
}