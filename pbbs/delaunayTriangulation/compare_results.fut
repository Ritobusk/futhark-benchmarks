-- > comparePoints ($loaddata "t10") ($loaddata "t20")

def comparePoints [n] (g1 : [n][2]i64) (g2: [n][2]i64) =
    let flags = map3 (\ a b r -> 
            if  a[0]== b[0] && a[1] ==b[1] then
                [-1i64, -1i64, -1i64, -1i64, -1i64]
            else
                [a[0], a[1],b[0], b[1], r]
        ) g1 g2 (indices g1)
    let diffs = filter (\x -> x[2] > -1) flags
    let t1 = trace flags
    in diffs
