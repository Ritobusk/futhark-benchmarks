--Shifts points' center towards (0,0)
def shiftPoints [n][d]
         (points : [n][d]f32) : [n][d]f32  =
     let mass_acc       = replicate (d) 0.0f32
     let mass_sum       = reduce (\acc point -> map2 (+) acc point ) mass_acc points
     let mass           = map (\i -> i / f32.i64 n) mass_sum
     let shifted_Points = map (\point -> map2 (\pv mv -> pv - mv) point mass) points
     in shifted_Points

def dist (p : (f32, f32)) (q: (f32, f32)) =
    f32.sqrt ((p.0 - q.0)**2 + (p.1 - q.1)**2)

def stencilK (k : i64) =
    [(k,0), (k,k),(k,-k),(0,k),(0,-k),(-k,0),(-k,k),(-k,-k)]

def log2Int (n : i64) : i64 =
  let (_, res) =
    loop (n, r) = (n, 0)
      while n > 1 do
        (n >> 1, r+1)
  in res 

type QorA = #left     | #right   | #up       | #down      | #center |
           #leftUp   | #rightUp | #leftDown | #rightDown

def findQuodrantOrAxis (p : (i64, i64)) (q : (i64, i64)) (grid_size : i64) : QorA =
    -- The input coordinates of a point will be (y,x) since its the row/col from the grid
    let (dir_x, dir_y) = (q.1 - p.1 ,(grid_size - q.0) - (grid_size -p.0))
    in if dir_x == 0 || dir_y == 0 then
        if      dir_x == 0 && dir_y > 0 then #up
        else if dir_x == 0 && dir_y < 0 then #down
        else if dir_x < 0 && dir_y == 0 then #left
        else if dir_x > 0 && dir_y == 0 then #right
        else #center
    else 
        if      dir_x > 0 && dir_y > 0 then #rightUp
        else if dir_x < 0 && dir_y > 0 then #leftUp
        else if dir_x < 0 && dir_y < 0 then #leftDown
        else #rightDown

def checkQuadrant (p : i32) (q1 : i32) (q2 : i32) (q3: i32) = p == q1 || p == q2 || p == q3

def classifyVertex (q1 : i32) (q2 : i32) (q3 : i32) (q4 : i32) : bool =
    if q1 != q3 && q2 != q4 then -- Diagonals are different
        if (q1 != q2 || q4 != q3) && (q1 != q4 || q2 != q3) then true -- Atleast 1 row or column has unique values
        else false
    else false

def classifyVertexI32 (q1 : i32) (q2 : i32) (q3 : i32) (q4 : i32) : i32 =
    if q1 != q3 && q2 != q4 then -- Diagonals are different
        if (q1 != q2 || q4 != q3) && (q1 != q4 || q2 != q3) then 1 -- Atleast 1 row or column has unique values
        else 0
    else 0

-- Returns triangles in clockwise order    
def classifyVertexAndTriangulation (q1 : i32) (q2 : i32) (q3 : i32) (q4 : i32) : ((i32, i32, i32),(i32, i32, i32)) = --: ((), ()) =
    if q1 != q3 && q2 != q4 then -- Diagonals are different
        if q1 != q2 && q1 != q4 && q2 != q3 && q4 != q3   then ((q2, q1, q4), (q4, q3, q2)) -- All are unique
        else if (q1 != q2 || q4 != q3) && (q1 != q4 || q2 != q3) then -- Atleast 1 row or column has unique values
            if q1 == q2 || q1 == q4 then ((q2, q4, q3), (-1, -1, -1))
            else if q2 == q1 || q2 == q3 then ((q1, q4, q3), (-1, -1, -1))
            else if q3 == q2 || q3 == q4 then ((q2, q1, q4), (-1, -1, -1))
            else ((q3, q2, q1), (-1, -1, -1))
        else ((-1, -1, -1), (-1, -1, -1))
    else ((-1, -1, -1), (-1, -1, -1))

def gridToGray [m] (grid : [m][m]i32) (max_idx : i64) : [m][m]f32  =
     tabulate_2d m m (\i j -> (f32.i32 grid[i][j]) / (f32.i64 max_idx))

def colours [m] (grid : [m][m]i32) : [m][m]u32 =
    let f (x) =
        let x = x**2
        in
            (u32.i32 (x*3) & 0xFF) << 16 |
            (u32.i32 (x+127) & 0xFF) << 8 |
            (u32.f64 (f64.cos (f64.i32 x) -f64.sin 3) & 0xFF)
    in map (map f) (grid)
