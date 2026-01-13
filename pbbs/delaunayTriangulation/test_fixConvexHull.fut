import "util"
import "delaunay_triangulation_VD"


def fixConvexHull2 [grid_size] [n] (grid : [grid_size][grid_size]i32) (points : [n][2]i64) = -- [](i32,i32,i32)
    let edge = map (\i -> 
        if i < grid_size then grid[0][i]
        else if i < 2 * grid_size then grid[(i % grid_size)][grid_size -1]
        else if i < 3 * grid_size then grid[grid_size -1][grid_size - 1 - (i % grid_size)]
        else grid[grid_size - 1 - (i% grid_size)][0]

        ) (iota (4 * grid_size)) 
    let edge = pack_points_i32 edge 
    
    let (_, triangles) = 
        loop (stack, triangles) = ([edge[0], edge[1]], [])
            for x in edge[2:] do 
                let stack = stack ++ [x]
                in if isClockwise points[stack[0]] points[stack[1]] points[stack[2]] grid_size then
                    (stack[1:], triangles)
                else
                    ([stack[0], stack[2]], triangles ++ [(stack[0], stack[2], stack[1])])
                
    in triangles
                

-- ==
-- entry: main
-- compiled random input {       [1000][2]f64 } 
-- compiled random input {    [1000000][2]f64 } 
def main [n]
    (points : [n][2]f64)  =
    let grid_size = 1024 

    let (grid, _, scaled_points) = movePointsToGrid points grid_size

    let grid' = voronoiDiagram  grid  scaled_points
    let grid''  = removeIslands grid' scaled_points
    let voronoi_diagram = tabulate_2d grid_size grid_size (\i j -> grid''[i][j].1) 

    let t1 = (fixConvexHull2 voronoi_diagram scaled_points) 
    let t1_len = length t1
    let t1 = sized (t1_len) t1
    let t2 = (fixConvexHull voronoi_diagram scaled_points) 
    let t2_len = length t2
    let t2 = sized (t1_len) t2
    
    
    let flags = map2 (\x y -> x.0 == y.0 && x.1 == y.1 && x.2 == y.2) t1 t2

    in reduce (&&) true flags && (t2_len == t1_len)
