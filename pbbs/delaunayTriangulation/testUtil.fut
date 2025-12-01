import "util"

def testClassifyVertex =
    let a = classifyVertex 0 1 2 3
    let b = classifyVertex 0 0 2 3
    let c = classifyVertex 0 1 1 3
    let d = classifyVertex 0 1 2 0
    let e = classifyVertex 0 1 2 2
    let satisfies = reduce (&&) true [a,b,c,d,e] 

    let f =  classifyVertex 0 1 2 1
    let g =  classifyVertex 0 1 0 3
    let h =  classifyVertex 1 1 0 0
    let i =  classifyVertex 0 1 1 0
    let j =  classifyVertex 0 1 1 1
    let k =  classifyVertex 1 1 1 0
    let not_satisfies = ! reduce (||) false [f,g,h,i,j, k] 
    in (satisfies, not_satisfies)

def main =

    let (t41, t42) = testClassifyVertex
    let t41 = trace t41
    let t42 = trace t42
    in 1
