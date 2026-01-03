--Shifts points' center towards (0,0)
def shiftPoints [n][d]
         (points : [n][d]f32) : [n][d]f32  =
     let mass_acc       = replicate (d) 0.0f32
     let mass_sum       = reduce (\acc point -> map2 (+) acc point ) mass_acc points
     let mass           = map (\i -> i / f32.i64 n) mass_sum
     let shifted_Points = map (\point -> map2 (\pv mv -> pv - mv) point mass) points
     in shifted_Points


def neq lte x y = if x `lte` y then !(y `lte` x) else true

def pack [n] lte (xs : [n](i64, i64)) =
  let (used, unused) =zip3 (indices xs) xs (rotate (-1) xs)
  |> partition (\(i,x,y) -> i == 0 || neq lte x.0 y.0) 
  in (map (.1) used, map (.1) unused)

def pack_points = pack (i64.<=)     


def dist (p : (f32, f32)) (q: (f32, f32)) : f32 =
    f32.sqrt ((p.0 - q.0)**2 + (p.1 - q.1)**2)

def disti64 (p : (i64, i64)) (q: (i64, i64)) : f32 =
    let p = (f32.i64 p.0, f32.i64 p.1)
    let q = (f32.i64 q.0, f32.i64 q.1)
    in f32.sqrt ((p.0 - q.0)**2 + (p.1 - q.1)**2)

def normalizei64 (p : (i64,i64)) : (f32, f32) =
    let p = (f32.i64 p.0, f32.i64 p.1)
    let mag = f32.sqrt (p.0**2 + p.1**2)
    in (p.0/mag, p.1/mag)

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
def classifyVertexAndTriangulation (q1 : i32) (q2 : i32) (q3 : i32) (q4 : i32) : ((i32, i32, i32),(i32, i32, i32)) = 
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
    
def triangleGrid [m] [k] [n] (grid : [m][m]i32) (voronoi_diagram : [m][m]i32) (triangles : [k](i32, i32, i32)) (points : [n][2]i64) = -- : [m][m]u32 =
    let colour_grid = colours voronoi_diagram
    let line_mask = loop (g) = (replicate m (replicate m 0)) for t in triangles do
            let (p1, p2, p3) = (points[t.0], points[t.1], points[t.2])
            let dir1 = normalizei64 (p2[0] - p1[0], p2[1] - p1[1])
            let dir2 = normalizei64 (p3[0] - p2[0], p3[1] - p2[1])
            let dir3 = normalizei64 (p1[0] - p3[0], p1[1] - p3[1])

            let dist1 = disti64 (p1[0], p1[1]) (p2[0], p2[1])
            let dist2 = disti64 (p2[0], p2[1]) (p3[0], p3[1])
            let dist3 = disti64 (p3[0], p3[1]) (p1[0], p1[1])

            let l1 = map (\i -> 
                let x = ((i64.f32 ((f32.i64 i) * dir1.0)) + p1[0])
                let y = ((i64.f32 ((f32.i64 i) * dir1.1)) + p1[1])
                in (y, x)
                ) <| iota (i64.f32 dist1) 
            let l2 = map (\i -> 
                let x = ((i64.f32 ((f32.i64 i) * dir2.0)) + p2[0])
                let y = ((i64.f32 ((f32.i64 i) * dir2.1)) + p2[1])
                in (y, x)
                ) <| iota (i64.f32 dist2) 
            let l3 = map (\i -> 
                let x = ((i64.f32 ((f32.i64 i) * dir3.0)) + p3[0])
                let y = ((i64.f32 ((f32.i64 i) * dir3.1)) + p3[1])
                in (y, x)
                ) <| iota (i64.f32 dist3) 

            let lines = l1 ++ l2 ++ l3
            in loop g' = g for x in lines do
                g' with [x.0,x.1] = 1
    in 
        tabulate_2d m m
            ( \r c ->
                if grid[r][c] >= 0 then 0u32
                else
                    if line_mask[r][c] > 0 then
                        u32.highest 
                    else
                        colour_grid[r][c]
            )
    
