use std::collections::{BTreeMap, BTreeSet};

use cgmath::prelude::*;
use cgmath::{point3, BaseFloat, Deg, Point3, Vector3};


type OrientedEdge = (usize, usize);


pub fn opposite(e: &OrientedEdge) -> OrientedEdge {
    (e.1, e.0)
}


#[derive(Debug)]
pub struct Mesh<T> {
    vertices: Vec<T>,
    at_vertex: Vec<OrientedEdge>,
    on_face: Vec<OrientedEdge>,
    on_boundary_component: Vec<OrientedEdge>,
    to_face: BTreeMap<OrientedEdge, usize>,
    next_on_face: BTreeMap<OrientedEdge, OrientedEdge>
}


impl<T> Mesh<T> {
    pub fn empty() -> Mesh<T> {
        Mesh {
            vertices: vec![],
            at_vertex: vec![],
            on_face: vec![],
            on_boundary_component: vec![],
            to_face: BTreeMap::new(),
            next_on_face: BTreeMap::new(),
        }
    }

    fn from_oriented_faces_unchecked(
        vertices: Vec<T>,
        face_lists: Vec<Vec<usize>>
    ) -> Mesh<T> {
        let oriented_edges_lists: Vec<_> = face_lists.iter()
            .map(cyclic_pairs).collect();

        let oriented_edges: Vec<_> = oriented_edges_lists.iter()
            .flatten().cloned().collect();

        let oriented_edge_set: BTreeSet<_> = oriented_edges.iter()
            .cloned().collect();

        let boundary_edges: Vec<_> = oriented_edges.iter()
            .filter(|e| !oriented_edge_set.contains(&opposite(e)))
            .cloned().collect();

        let boundary_lists = boundary_cycles(boundary_edges);

        let at_vertex: Vec<_> = oriented_edges.iter()
            .map(|&(v, w)| (v, (v, w)))
            .collect::<BTreeMap<_, _>>().values().cloned()
            .collect();

        let to_face: BTreeMap<_, _> = oriented_edges_lists.iter()
            .enumerate()
            .flat_map(|(i, f)| f.iter().map(move |&e| (e, i)))
            .collect();

        let along_face: Vec<_> = to_face.iter()
            .map(|(&e, &f)| (f, e))
            .collect::<BTreeMap<_, _>>().iter()
            .map(|(_, &e)| e)
            .collect();

        let to_boundary_component: BTreeMap<_, _> = boundary_lists.iter()
            .enumerate()
            .flat_map(|(i, b)| b.iter().map(move |&e| (e, i)))
            .collect();

        let on_boundary_component: Vec<_> = to_boundary_component.iter()
            .map(|(&e, &b)| (b, e))
            .collect::<BTreeMap<_, _>>().iter()
            .map(|(_, &e)| e)
            .collect();

        let next_on_face: BTreeMap<_, _> = oriented_edges_lists.iter()
            .chain(boundary_lists.iter())
            .flat_map(cyclic_pairs)
            .collect();

        Mesh {
            vertices,
            at_vertex,
            on_face: along_face,
            on_boundary_component,
            to_face,
            next_on_face,
        }
    }

    pub fn from_oriented_faces<
        I1: IntoIterator<Item=T>,
        I2: IntoIterator<Item=usize>,
        I3: IntoIterator<Item=I2>
    >(
        vertices: I1,
        face_lists_in: I3
    ) -> Result<Mesh<T>, String> {
        let vertices: Vec<_> = vertices.into_iter().collect();
        let face_lists: Vec<Vec<_>> = face_lists_in.into_iter()
            .map(|f| f.into_iter().collect())
            .collect();

        let defined_vertices: BTreeSet<_> = (0..vertices.len()).collect();

        let seen_vertices: BTreeSet<_> = face_lists.iter()
            .filter(|f| f.len() > 0)
            .flatten().cloned().collect();

        let oriented_edges: Vec<_> = face_lists.iter()
            .filter(|f| f.len() > 0)
            .map(cyclic_pairs)
            .flatten().collect();

        let oriented_edge_set: BTreeSet<_> = oriented_edges.iter()
            .cloned().collect();

        let boundary_vertices: Vec<_> = oriented_edges.iter()
            .filter(|e| !oriented_edge_set.contains(&opposite(e)))
            .map(|&(v, _)| v)
            .collect();

        if seen_vertices.iter().any(|v| !defined_vertices.contains(v)) {
            Err("an undefined vertex appears in a face".to_string())
        } else if defined_vertices.iter().any(|v| !seen_vertices.contains(v)) {
            Err("some vertex does not appear in any faces".to_string())
        } else if face_lists.iter().any(|f| f.len() < 2) {
            Err("some face has fewer than two vertices".to_string())
        } else if face_lists.iter().any(|f| !all_unique(f)) {
            Err("a vertex appears more than once in the same face".to_string())
        } else if !all_unique(boundary_vertices) {
            Err("a vertex appears more then once in a boundary".to_string())
        } else if !all_unique(oriented_edges) {
            Err("an oriented edge appears more than once".to_string())
        } else {
            Ok(Mesh::from_oriented_faces_unchecked(vertices, face_lists))
        }
    }

    pub fn vertices<'a>(&'a self) -> &'a Vec<T> {
        &self.vertices
    }

    pub fn previous_at_vertex(&self, e: OrientedEdge) -> OrientedEdge {
        self.next_on_face.get(&opposite(&e)).unwrap().clone()
    }

    pub fn next_on_face(&self, e: OrientedEdge) -> OrientedEdge {
        self.next_on_face.get(&e).unwrap().clone()
    }

    pub fn vertices_in_face(&self, start: OrientedEdge) -> Vec<usize> {
        canonical_circular(
            trace_cycle(start, |e| Some(self.next_on_face(e)))
                .iter()
                .map(|&(v, _)| v)
                .collect()
        )
    }

    pub fn face_indices(&self) -> Vec<Vec<usize>> {
        self.on_face.iter()
            .map(|&e| self.vertices_in_face(e))
            .collect()
    }

    pub fn boundary_indices(&self) -> Vec<Vec<usize>> {
        self.on_boundary_component.iter()
            .map(|&e| self.vertices_in_face(e))
            .collect()
    }

    pub fn edge_indices(&self) -> Vec<(usize, usize)> {
        self.next_on_face.keys()
            .map(|&(u, v)| (u.min(v), u.max(v)))
            .collect::<BTreeSet<_>>().into_iter()
            .collect()
    }

    pub fn boundary_edge_indices(&self) -> Vec<(usize, usize)> {
        self.boundary_indices().iter()
            .flat_map(cyclic_pairs)
            .collect()
    }

    pub fn vertex_neighbors(&self, start: OrientedEdge) -> Vec<usize> {
        canonical_circular(
            trace_cycle(start, |e| Some(self.previous_at_vertex(e)))
                .iter()
                .map(|&(_, w)| w)
                .rev()
                .collect()
        )
    }

    pub fn neighbor_indices(&self) -> Vec<Vec<usize>> {
        self.at_vertex.iter()
            .map(|&e| self.vertex_neighbors(e))
            .collect()
    }

    pub fn map_vertices<S>(&self, f: impl Fn(&T) -> S)
        -> Result<Mesh<S>, String>
    {
        Mesh::from_oriented_faces(
            self.vertices.iter().map(f),
            self.face_indices(),
        )
    }
}


impl<T: Clone> Mesh<T> {
    pub fn triangulate(&self) -> Result<Self, String> {
        Self::from_oriented_faces(
            self.vertices.clone(),
            self.face_indices().iter().flat_map(triangulate)
        )
    }

    pub fn quadrangulate(&self, compose: impl Fn(&Vec<T>) -> T) -> Self {
        let vertices = self.vertices();
        let edges = self.edge_indices();
        let nr_vertices = vertices.len();
        let nr_edges = edges.len();

        let mid_point_index: BTreeMap<_, _> = edges.iter().enumerate()
            .map(|(i, e)| (e, i + nr_vertices))
            .flat_map(|(&(u, v), i)| [((u, v), i), ((v, u), i)])
            .collect();

        let verts_out: Vec<_> = (0..nr_vertices).map(|i| vec![i])
            .chain(edges.iter().map(|&(u, v)| vec![u, v]))
            .chain(self.face_indices())
            .map(|is| compose(
                &is.iter().map(|&i| vertices[i].clone()).collect()
            ))
            .collect();

        let mut faces_out = vec![];
        for (&(u, v), &i) in self.to_face.iter() {
            let w = self.next_on_face((u, v)).1;

            faces_out.push(vec![
                v,
                mid_point_index[&(v, w)],
                i + nr_vertices + nr_edges,
                mid_point_index[&(u, v)]
            ]);
        }

        Self::from_oriented_faces_unchecked(verts_out, faces_out)
    }
}


impl<S: BaseFloat> Mesh<Point3<S>> {
    pub fn subd(&self, fix_boundary: bool) -> Self {
        let nr_verts = self.vertices().len();
        let nr_edges = self.edge_indices().len();
        let sub_mesh = self.quadrangulate(|vs| Point3::centroid(vs));
        let neighbor_indices = sub_mesh.neighbor_indices();
        let vertices_tmp = sub_mesh.vertices();

        let mut vertices_out = sub_mesh.vertices().clone();

        let bnd: BTreeSet<_> = sub_mesh.boundary_indices().iter()
            .flatten().cloned().collect();

        // adjust edge-point positions for those not on the boundary
        for k in nr_verts..(nr_verts + nr_edges) {
            if !bnd.contains(&k) {
                let nbs = neighbor_indices[k].clone();
                vertices_out[k] = indexed_centroid(&vertices_tmp, nbs);
            }
        }

        // adjust the positions of the original vertices
        for k in 0..nr_verts {
            let pos_in = vertices_out[k];
            let nbs = neighbor_indices[k].clone();

            if !bnd.contains(&k) {
                let c1 = indexed_centroid(&vertices_tmp, nbs.clone());
                let c2 = indexed_centroid(&vertices_out, nbs.clone());

                vertices_out[k] = Point3::centroid(
                    &vec![pos_in; nbs.len() - 3].iter()
                        .chain([c1, c2, c2].iter())
                        .cloned().collect::<Vec<_>>()
                );
            } else if !fix_boundary {
                let nbs_bnd = nbs.iter().filter(|v| bnd.contains(v)).cloned();
                let c = indexed_centroid(&vertices_tmp, nbs_bnd);
                vertices_out[k] = Point3::centroid(&[pos_in, c]);
            }
        }

        let faces = sub_mesh.face_indices();
        Mesh::from_oriented_faces_unchecked(vertices_out, faces)
    }

    pub fn tightened(&self, fix_boundary: bool) -> Self {
        let pos = self.vertices();
        let face_indices = self.face_indices();
        let n = pos.len();

        let mut px = vec![0.0; n];
        let mut py = vec![0.0; n];
        let mut pz = vec![0.0; n];
        for v in 0..n {
            px[v] = pos[v].x.to_f64().unwrap();
            py[v] = pos[v].y.to_f64().unwrap();
            pz[v] = pos[v].z.to_f64().unwrap();
        }

        let mut fixed = vec![false; n];
        if fix_boundary {
            for bnd in self.boundary_indices() {
                for v in bnd {
                    fixed[v] = true;
                }
            }
        }

        let mut gx = vec![0.0; n];
        let mut gy = vec![0.0; n];
        let mut gz = vec![0.0; n];

        let f1 = 0.1;
        let f2 = 0.1;

        for _ in 0..100 {
            for v in 0..n {
                gx[v] = 0.0;
                gy[v] = 0.0;
                gz[v] = 0.0;
            }

            for f in &face_indices {
                let m = f.len();

                for k in 0..m {
                    let u = f[k];
                    let v = f[(k + 1) % m];
                    let w = f[(k + 2) % m];

                    let ax = px[u] - px[v];
                    let ay = py[u] - py[v];
                    let az = pz[u] - pz[v];

                    let bx = px[w] - px[v];
                    let by = py[w] - py[v];
                    let bz = pz[w] - pz[v];

                    let cx = px[w] - px[u];
                    let cy = py[w] - py[u];
                    let cz = pz[w] - pz[u];

                    let nx = by * az - bz * ay;
                    let ny = bz * ax - bx * az;
                    let nz = bx * ay - by * ax;

                    let nl = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-6);

                    gx[v] += (ny * cz - nz * cy) / nl + f1 * (ax + bx);
                    gy[v] += (nz * cx - nx * cz) / nl + f1 * (ay + by);
                    gz[v] += (nx * cy - ny * cx) / nl + f1 * (az + bz);
                }
            }

            for v in 0..n {
                if !fixed[v] {
                    px[v] += f2 * gx[v];
                    py[v] += f2 * gy[v];
                    pz[v] += f2 * gz[v];
                }
            }
        }

        let mut pos_out = vec![];
        for v in 0..n {
            pos_out.push(point3(
                S::from(px[v]).unwrap(),
                S::from(py[v]).unwrap(),
                S::from(pz[v]).unwrap(),
            ))
        }

        Mesh::from_oriented_faces_unchecked(pos_out, self.face_indices())
    }

    pub fn inset_corner(&self, bnd_edge: OrientedEdge, wd: S) -> Point3<S> {
        let pos = self.vertices();
        let mut e = bnd_edge;
        let mut ends = vec![];

        loop {
            ends.push(pos[e.1]);
            e = self.previous_at_vertex(e);
            if e == bnd_edge {
                break;
            }
        }

        let corner = pos[bnd_edge.0];
        let left = ends[0];
        let right = ends[ends.len() - 1];

        let points = if ends.len() > 2 {
            &ends[1..(ends.len() - 1)]
        } else {
            &self.vertices_in_face(self.previous_at_vertex(bnd_edge)).iter()
                .map(|&v| pos[v])
                .collect::<Vec<_>>()[..]
        };

        let center = Point3::centroid(points);

        inset_point(corner, wd, left, right, center)
    }

    pub fn elevate_corner(&self, bnd_edge: OrientedEdge, ht: S) -> Point3<S> {
        let pos = self.vertices();
        let mut ea = bnd_edge;
        let mut n = Vector3::zero();

        loop {
            let eb = self.previous_at_vertex(ea);
            if eb == bnd_edge {
                break;
            }
            n += (pos[eb.1] - pos[eb.0]).cross(pos[ea.1] - pos[ea.0]);
            ea = eb;
        }

        pos[bnd_edge.0] + n.normalize() * ht
    }

    pub fn revised_boundaries(&self, f: impl Fn(OrientedEdge) -> Point3<S>)
        -> Mesh<Point3<S>>
    {
        let mut pos = self.vertices().clone();

        for &e in &self.on_boundary_component {
            let mods = self.cycle_revisions(&f, e);
            for (v, p) in mods {
                pos[v] = p;
            }
        }

        Mesh::from_oriented_faces(pos, self.face_indices()).unwrap()
    }

    pub fn boundary_strips(&self, f: impl Fn(OrientedEdge) -> Point3<S>)
        -> Mesh<Point3<S>>
    {
        let mut pos = vec![];
        let mut faces = vec![];

        for &e in &self.on_boundary_component {
            let mods = self.cycle_revisions(&f, e);
            let n = pos.len();
            let m = mods.len();

            pos.extend(mods.iter().map(|&(v, _)| self.vertices()[v]));
            pos.extend(mods.iter().map(|&(_, p)| p));

            for i in 0..m {
                let j = (i + 1) % m;
                faces.push(vec![n + m + i, n + m + j, n + j, n + i]);
            }
        }

        Mesh::from_oriented_faces(pos, faces).unwrap()
    }

    fn cycle_revisions(
        &self,
        f: impl Fn(OrientedEdge) -> Point3<S>,
        start: OrientedEdge
    )
        -> Vec<(usize, Point3<S>)>
    {
        let mut result = vec![];
        let mut e = start;

        loop {
            result.push((e.0, f(e)));
            e = self.next_on_face(e);
            if e == start {
                break;
            }
        }

        result
    }
}


fn boundary_cycles(boundary_edges: Vec<OrientedEdge>)
    -> Vec<Vec<(usize, usize)>>
{
    let items: Vec<_> = boundary_edges.iter().map(|&(v, _)| v).collect();
    let advance: BTreeMap<_, _> = boundary_edges.iter().map(opposite).collect();

    let mut seen: BTreeSet<usize> = BTreeSet::new();
    let mut result = vec![];

    for v in items {
        if !seen.contains(&v) {
            let cycle = trace_cycle(v, |v| advance.get(&v).copied());
            seen.extend(&cycle);
            result.push(cyclic_pairs(&cycle));
        }
    }

    result
}


fn trace_cycle<T: Copy + PartialEq>(v: T, advance: impl Fn(T) -> Option<T>)
    -> Vec<T>
{
    let mut cycle = vec![];
    let mut w = v;

    loop {
        if let Some(u) = advance(w) {
            cycle.push(w);
            w = u;
            if w == v {
                break;
            }
        } else {
            cycle.clear();
            break;
        }
    }

    cycle
}


fn triangulate<T: Copy>(corners: &Vec<T>) -> Vec<Vec<T>> {
    (1..(corners.len() - 1))
        .map(|i| vec![corners[0], corners[i], corners[i + 1]])
        .collect()
}


fn cyclic_pairs<T: Copy>(indices: &Vec<T>) -> Vec<(T, T)> {
    let mut result = vec![];
    for i in 0..(indices.len() - 1) {
        result.push((indices[i], indices[i + 1]))
    }
    result.push((indices[indices.len() - 1], indices[0]));
    result
}


fn canonical_circular<T: Copy + Ord>(list: Vec<T>) -> Vec<T> {
    if list.is_empty() {
        vec![]
    } else {
        (0..list.len()).map(|k|
            list.iter().skip(k).chain(list.iter().take(k)).cloned().collect()
        ).min().unwrap()
    }
}


fn all_unique<T: Ord, I: IntoIterator<Item=T>>(items: I) -> bool {
    let mut seen: BTreeSet<T> = BTreeSet::new();

    for x in items.into_iter() {
        if seen.contains(&x) {
            return false;
        }
        seen.insert(x);
    }

    true
}


fn indexed_centroid<S: BaseFloat, I: IntoIterator<Item=usize>>(
    positions: &Vec<Point3<S>>,
    indices: I
) -> Point3<S>
{
    Point3::centroid(
        &indices.into_iter()
            .map(|v| positions[v])
            .collect::<Vec<_>>()
    )
}


pub fn inset_point<S: BaseFloat>(
    corner: Point3<S>,
    wd: S,
    left: Point3<S>,
    right: Point3<S>,
    center: Point3<S>,
) -> Point3<S>
{
    let norm = |v: Vector3<S>| v.dot(v).sqrt();
    let normx = |v: Vector3<S>| norm(v).max(S::from(1e-6).unwrap());
    let normalized = |v: Vector3<S>| v / normx(v);
    let project_along = |v: Vector3<S>, n: Vector3<S>| v - n * v.dot(n);

    let lft = normalized(left - corner);
    let rgt = normalized(right - corner);
    let dia = lft + rgt;

    if norm(dia).to_f64().unwrap() < 0.01 {
        let s = normalized(center - corner);
        let t = project_along(s, lft);
        return corner + s * (wd / normx(t));
    } else if norm(lft.cross(rgt)).to_f64().unwrap() < 0.01 {
        return corner + normalized(dia) * wd;
    } else {
        let len = wd * norm(dia) / normx(project_along(dia, lft));
        let s = normalized(dia.cross(lft.cross(rgt)));
        let t = project_along(center - corner, s);
        let f = len / normx(t);
        return corner + t * f;
    }
}


pub fn cylinder<S: BaseFloat>(
    start: Point3<S>, end: Point3<S>, radius: S, nr_segments: usize
)
    -> Mesh<Point3<S>>
{
    let dir = (end - start).normalize();

    let t: Vector3<S> =
        if Vector3::unit_x().dot(dir).abs().to_f64().unwrap() < 0.71 {
            Vector3::unit_x()
        } else {
            Vector3::unit_y()
        };

    let a = (t - dir * t.dot(dir)).normalize();
    let b = dir.cross(a).normalize();

    assert!((a.cross(b).dot(dir).to_f64().unwrap() - 1.0).abs() < 1e-12);

    let mut pos = vec![];
    let mut faces = vec![];

    for p in [start, end] {
        for i in 0..nr_segments {
            let alpha = Deg(S::from(360.0 / 24.0 * i as f64).unwrap());
            pos.push(p + (a * alpha.cos() + b * alpha.sin()) * radius);
        }
    }

    for i in 0..nr_segments {
        let j = (i + 1) % nr_segments;
        faces.push(vec![nr_segments + i, nr_segments + j, j, i]);
    }

    Mesh::from_oriented_faces(pos, faces).unwrap()
}


fn spherified<S: BaseFloat>(mesh: &Mesh<Point3<S>>) -> Mesh<Point3<S>> {
    Mesh::from_oriented_faces_unchecked(
        mesh.vertices().iter()
            .map(|p| Point3::origin() + (p - Point3::origin()).normalize())
            .collect(),
        mesh.face_indices()
    )
}


pub fn unit_sphere<S: BaseFloat>(nr_divisions: usize)
    -> Mesh<Point3<S>>
{
    let one = S::one();

    let mut mesh = Mesh::from_oriented_faces(
        [
            point3( one,  one,  one), // 0
            point3(-one, -one,  one), // 1
            point3( one, -one, -one), // 2
            point3(-one,  one, -one), // 3
        ],
        [
            [0, 1, 2],
            [0, 2, 3],
            [0, 3, 1],
            [3, 2, 1],
        ]
    )
        .unwrap();

    mesh = spherified(&mesh);

    for _ in 0..nr_divisions {
        mesh = spherified(&mesh.subd(true))
    }

    mesh
}


pub enum ItemType {
    Vertex,
    Edge,
    Face,
}


pub fn decompose_mesh<S: BaseFloat>(mesh: Mesh<Point3<S>>)
    -> Vec<(Mesh<Point3<S>>, ItemType)>
{
    let edge_radius = S::from(0.049).unwrap();
    let edge_lift = S::from(0.05).unwrap();

    let mut result = vec![];

    for (u, v) in mesh.edge_indices() {
        let start = mesh.vertices()[u];
        let end = mesh.vertices()[v];
        let edge = cylinder(start, end, edge_radius, 24);

        result.push((edge, ItemType::Edge));
    }

    let s = unit_sphere(4);

    for center in mesh.vertices() {
        let sphere = Mesh::from_oriented_faces_unchecked(
            s.vertices().iter()
                .map(|p| center + (p - Point3::origin()) * edge_radius)
                .collect(),
            s.face_indices()
        );
        result.push((sphere, ItemType::Vertex));
    }

    for f in mesh.face_indices() {
        let mut face_mesh = Mesh::from_oriented_faces(
            f.iter().map(|&v| mesh.vertices()[v]),
            [(0..f.len())]
        ).unwrap();

        for _ in 0..4 {
            face_mesh = face_mesh.subd(true).tightened(true);
        }

        let elevate = |e| face_mesh.elevate_corner(e, edge_lift);
        let elevated = face_mesh.revised_boundaries(elevate).tightened(true);

        result.push((elevated, ItemType::Face));
    }

    result
}


#[cfg(test)]
mod test {
    use cgmath::{point3, EuclideanSpace, Point3};

    use super::*;

    fn octa_vert_names() -> [String; 6] {
        [
            "front".to_string(),
            "right".to_string(),
            "top".to_string(),
            "back".to_string(),
            "left".to_string(),
            "bottom".to_string(),
        ]
    }

    fn octa_vert_pos() -> [Point3<f32>; 6] {
        [
            point3( 1.0,  0.0,  0.0),
            point3( 0.0,  1.0,  0.0),
            point3( 0.0,  0.0,  1.0),
            point3(-1.0,  0.0,  0.0),
            point3( 0.0, -1.0,  0.0),
            point3( 0.0,  0.0, -1.0),
        ]
    }

    fn octa_faces() -> [Vec<usize>; 8] {
        [
            vec![ 0, 1, 2 ],
            vec![ 1, 0, 5 ],
            vec![ 2, 1, 3 ],
            vec![ 0, 2, 4 ],
            vec![ 3, 5, 4 ],
            vec![ 5, 3, 1 ],
            vec![ 4, 5, 0 ],
            vec![ 3, 4, 2 ],
        ]
    }

    #[test]
    fn test_cyclic_pairs() {
        assert_eq!(cyclic_pairs(&vec![1, 2, 3]), vec![(1, 2), (2, 3), (3, 1)]);
    }

    #[test]
    fn test_empty() {
        let mesh = Mesh::<i32>::empty();

        assert_eq!(mesh.vertices(), &[]);
        assert_eq!(mesh.edge_indices(), []);
        assert_eq!(mesh.face_indices(), [[]; 0]);
        assert_eq!(mesh.boundary_indices(), [[]; 0]);
        assert_eq!(mesh.neighbor_indices(), [[]; 0]);
    }

    #[test]
    fn test_without_boundary() {
        let octa = Mesh::from_oriented_faces(
            octa_vert_names(), octa_faces()
        ).unwrap();

        assert_eq!(
            octa.edge_indices(),
            [
                ( 0, 1 ),
                ( 0, 2 ),
                ( 0, 4 ),
                ( 0, 5 ),
                ( 1, 2 ),
                ( 1, 3 ),
                ( 1, 5 ),
                ( 2, 3 ),
                ( 2, 4 ),
                ( 3, 4 ),
                ( 3, 5 ),
                ( 4, 5 ),
            ]
        );

        assert_eq!(
            octa.face_indices().into_iter()
                .collect::<BTreeSet<_>>().into_iter()
                .collect::<Vec<_>>(),
            [
                [ 0, 1, 2 ],
                [ 0, 2, 4 ],
                [ 0, 4, 5 ],
                [ 0, 5, 1 ],
                [ 1, 3, 2 ],
                [ 1, 5, 3 ],
                [ 2, 3, 4 ],
                [ 3, 5, 4 ],
            ]
        );

        assert_eq!(octa.boundary_indices(), [[]; 0]);

        assert_eq!(
            octa.neighbor_indices(),
            [
                [ 1, 2, 4, 5 ],
                [ 0, 5, 3, 2 ],
                [ 0, 1, 3, 4 ],
                [ 1, 5, 4, 2 ],
                [ 0, 2, 3, 5 ],
                [ 0, 4, 3, 1 ],
            ]
        );
    }

    #[test]
    fn test_with_boundary() {
        let octa = Mesh::from_oriented_faces(
            octa_vert_names(),
            [
                [ 1, 0, 5 ],
                [ 2, 1, 3 ],
                [ 0, 2, 4 ],
                [ 5, 3, 1 ],
                [ 4, 5, 0 ],
                [ 3, 4, 2 ],
            ],
        ).unwrap();

        assert_eq!(
            octa.edge_indices(),
            [
                ( 0, 1 ),
                ( 0, 2 ),
                ( 0, 4 ),
                ( 0, 5 ),
                ( 1, 2 ),
                ( 1, 3 ),
                ( 1, 5 ),
                ( 2, 3 ),
                ( 2, 4 ),
                ( 3, 4 ),
                ( 3, 5 ),
                ( 4, 5 ),
            ]
        );

        assert_eq!(
            octa.face_indices().into_iter()
                .collect::<BTreeSet<_>>().into_iter()
                .collect::<Vec<_>>(),
            [
                [ 0, 2, 4 ],
                [ 0, 4, 5 ],
                [ 0, 5, 1 ],
                [ 1, 3, 2 ],
                [ 1, 5, 3 ],
                [ 2, 3, 4 ],
            ]
        );

        assert_eq!(octa.boundary_indices(), [[ 0, 1, 2 ], [ 3, 5, 4 ]]);

        assert_eq!(
            octa.neighbor_indices(),
            [
                [ 1, 2, 4, 5 ],
                [ 0, 5, 3, 2 ],
                [ 0, 1, 3, 4 ],
                [ 1, 5, 4, 2 ],
                [ 0, 2, 3, 5 ],
                [ 0, 4, 3, 1 ],
            ]
        );
    }

    #[test]
    fn test_undefined_vertex() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names()[1..].into_iter(),
                octa_faces()
            ).is_err()
        );
    }

    #[test]
    fn test_unreferenced_vertex() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names().into_iter().chain(["off".to_string()]),
                octa_faces()
            ).is_err()
        );
    }

    #[test]
    fn test_vertex_duplicate_in_face() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names(),
                [
                    vec![ 0, 1, 2, 0, 4, 5 ],
                    vec![ 1, 0, 5 ],
                    vec![ 2, 1, 3 ],
                    vec![ 0, 2, 4 ],
                    vec![ 3, 5, 4 ],
                    vec![ 5, 3, 1 ],
                    vec![ 3, 4, 2 ],
                ]
            ).is_err()
        );
    }

    #[test]
    fn test_vertex_duplicate_in_boundary() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names(),
                [
                    [ 1, 0, 5 ],
                    [ 2, 1, 3 ],
                    [ 0, 2, 4 ],
                    [ 3, 5, 4 ],
                    [ 5, 3, 1 ],
                    [ 3, 4, 2 ],
                ]
            ).is_err()
        );
    }

    #[test]
    fn test_empty_face() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names(),
                octa_faces().into_iter().chain([vec![]])
            ).is_err()
        );
    }

    #[test]
    fn test_one_gon() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names(),
                octa_faces().into_iter().chain([vec![0]])
            ).is_err()
        );
    }

    #[test]
    fn test_single_two_gon() {
        assert!(
            Mesh::from_oriented_faces(['a', 'b'], [[0, 1]]).is_ok()
        );
    }

    #[test]
    fn test_orientation_mismatch() {
        assert!(
            Mesh::from_oriented_faces(
                octa_vert_names(),
                octa_faces().into_iter().skip(1).chain([vec![0, 2, 1]])
            ).is_err()
        );
    }

    #[test]
    fn test_duplicate_edge() {
        assert!(
            Mesh::from_oriented_faces(
                ['a', 'b', 'c', 'd'],
                [
                    [ 2, 0, 1 ],
                    [ 0, 2, 3 ],
                    [ 2, 0, 3 ],
                    [ 0, 2, 1 ],
                ]
                ).is_err()
        );
    }

    #[test]
    fn test_vertices_method() {
        let mesh = Mesh::from_oriented_faces(
            octa_vert_names(),
            octa_faces()
        ).unwrap();

        assert_eq!(mesh.vertices(), &octa_vert_names());
    }

    #[test]
    fn test_triangulate() {
        assert_eq!(
            triangulate(&vec![0, 1, 2, 3]),
            vec![vec![0, 1, 2], vec![0, 2, 3]]
        );

        assert_eq!(
            triangulate(&vec![3, 2, 1, 0]),
            vec![vec![3, 2, 1], vec![3, 1, 0]]
        );
    }

    #[test]
    fn test_triangulate_mesh() {
        let mesh = Mesh::from_oriented_faces(
            ['a', 'b', 'c', 'd', 'e', 'f'],
            [
                vec![2, 1, 0],
                vec![3, 4, 5],
                vec![0, 1, 4, 3],
                vec![1, 2, 5, 4],
                vec![2, 0, 3, 5],
            ]
        ).unwrap().triangulate().unwrap();

        assert_eq!(mesh.vertices(), &['a', 'b', 'c', 'd', 'e', 'f']);
        assert_eq!(
            mesh.face_indices().into_iter()
                .collect::<BTreeSet<_>>().into_iter()
                .collect::<Vec<_>>(),
            [
                [0, 1, 4],
                [0, 2, 1],
                [0, 3, 5],
                [0, 4, 3],
                [0, 5, 2],
                [1, 2, 5],
                [1, 5, 4],
                [3, 4, 5],
            ]
        );
    }

    #[test]
    fn test_map_vertices() {
        let octa = Mesh::from_oriented_faces(
            octa_vert_names(), octa_faces()
        ).unwrap();

        let mapped = octa.map_vertices(|v| v[..2].to_uppercase()).unwrap();

        assert_eq!(mapped.vertices(), &["FR", "RI", "TO", "BA", "LE", "BO"]);
    }

    #[test]
    fn test_quadrangulate() {
        let octa = Mesh::from_oriented_faces(
            octa_vert_pos().iter().map(|p| p * 6.0).collect::<Vec<_>>(),
            octa_faces()
        ).unwrap();

        let sub = octa.quadrangulate(|vs| Point3::centroid(vs));

        assert_eq!(sub.vertices().len(), 26);
        assert_eq!(sub.edge_indices().len(), 48);
        assert_eq!(sub.face_indices().len(), 24);

        assert_eq!(
            sub.face_indices().iter()
                .map(Vec::len)
                .collect::<Vec<_>>(),
            &[4; 24]
        );

        assert_eq!(
            sub.face_indices().iter()
                .map(|f| f.iter().filter(|&&v| v < 6).count())
                .collect::<Vec<_>>(),
            &[1; 24]
        );

        assert_eq!(
            sub.face_indices().iter()
                .map(|f| f.iter().filter(|&&v| v < 18).count())
                .collect::<Vec<_>>(),
            &[3; 24]
        );

        assert_eq!(
            sub.neighbor_indices().iter().filter(|ns| ns.len() == 4).count(),
            18
        );

        assert_eq!(
            sub.neighbor_indices().iter().filter(|ns| ns.len() == 3).count(),
            8
        );

        for p in [
            point3( 0.0,  0.0,  6.0),
            point3( 0.0,  3.0, -3.0),
            point3( 0.0,  3.0,  3.0),
            point3( 0.0,  6.0,  0.0),
            point3( 2.0, -2.0, -2.0),
            point3( 2.0, -2.0,  2.0),
            point3( 2.0,  2.0, -2.0),
            point3( 2.0,  2.0,  2.0),
            point3( 3.0, -3.0,  0.0),
            point3( 3.0,  0.0, -3.0),
            point3( 3.0,  0.0,  3.0),
            point3( 3.0,  3.0,  0.0),
            point3( 6.0,  0.0,  0.0),
        ] {
            assert!(sub.vertices().contains(&p));
            assert!(sub.vertices().contains(&(p * -1.0)));
        }
    }

    #[test]
    fn test_subd() {
        let octa = Mesh::from_oriented_faces(
            octa_vert_pos().iter().map(|p| p * 12.0).collect::<Vec<_>>(),
            octa_faces()
        ).unwrap();

        let quad = octa.quadrangulate(|vs| Point3::centroid(vs));
        let subd = octa.subd(false);

        assert_eq!(quad.face_indices(), subd.face_indices());

        for (i, ns) in subd.neighbor_indices().iter().enumerate() {
            println!("{i}: {ns:?}");
        }
        println!();

        for (i, v) in subd.vertices().iter().enumerate() {
            println!("{i}: {v:?}");
        }

        for p in [
            point3( 0.0,  0.0,  7.0),
            point3( 0.0,  5.0, -5.0),
            point3( 0.0,  5.0,  5.0),
            point3( 0.0,  7.0,  0.0),
            point3( 4.0, -4.0, -4.0),
            point3( 4.0, -4.0,  4.0),
            point3( 4.0,  4.0, -4.0),
            point3( 4.0,  4.0,  4.0),
            point3( 5.0, -5.0,  0.0),
            point3( 5.0,  0.0, -5.0),
            point3( 5.0,  0.0,  5.0),
            point3( 5.0,  5.0,  0.0),
            point3( 7.0,  0.0,  0.0),
        ] {
            assert!(subd.vertices().contains(&p));
            assert!(subd.vertices().contains(&(p * -1.0)));
        }
    }

    #[test]
    fn test_subd_with_boundary_1() {
        let mesh = Mesh::from_oriented_faces(
            [
                point3( 0.0, 0.0, 0.0), // 0
                point3( 8.0, 0.0, 0.0), // 1
                point3(16.0, 0.0, 0.0), // 2
                point3( 0.0, 8.0, 0.0), // 3
                point3( 8.0, 8.0, 0.0), // 4
                point3(16.0, 8.0, 0.0), // 5
            ],
            [
                [0, 1, 4, 3],
                [1, 2, 5, 4],
            ]
        ).unwrap();

        let quad = mesh.quadrangulate(|vs| Point3::centroid(vs));
        let subd = mesh.subd(false);

        assert_eq!(quad.face_indices(), subd.face_indices());
        assert_eq!(subd.vertices().len(), 15);

        for p in [
            point3( 0.0, 4.0, 0.0),
            point3( 1.0, 1.0, 0.0),
            point3( 1.0, 7.0, 0.0),
            point3( 4.0, 0.0, 0.0),
            point3( 4.0, 4.0, 0.0),
            point3( 4.0, 8.0, 0.0),
            point3( 8.0, 0.0, 0.0),
            point3( 8.0, 4.0, 0.0),
            point3( 8.0, 8.0, 0.0),
            point3(12.0, 0.0, 0.0),
            point3(12.0, 4.0, 0.0),
            point3(12.0, 8.0, 0.0),
            point3(15.0, 1.0, 0.0),
            point3(15.0, 7.0, 0.0),
            point3(16.0, 4.0, 0.0),
        ] {
            assert!(subd.vertices().contains(&p));
        }
    }

    #[test]
    fn test_subd_with_fixed_boundary_1() {
        let mesh = Mesh::from_oriented_faces(
            [
                point3( 0.0, 0.0, 0.0), // 0
                point3( 8.0, 0.0, 0.0), // 1
                point3(16.0, 0.0, 0.0), // 2
                point3( 0.0, 8.0, 0.0), // 3
                point3( 8.0, 8.0, 0.0), // 4
                point3(16.0, 8.0, 0.0), // 5
            ],
            [
                [0, 1, 4, 3],
                [1, 2, 5, 4],
            ]
        ).unwrap();

        let quad = mesh.quadrangulate(|vs| Point3::centroid(vs));
        let subd = mesh.subd(true);

        assert_eq!(quad.face_indices(), subd.face_indices());
        assert_eq!(subd.vertices().len(), 15);

        for p in [
            point3( 0.0, 4.0, 0.0),
            point3( 0.0, 0.0, 0.0),
            point3( 0.0, 8.0, 0.0),
            point3( 4.0, 0.0, 0.0),
            point3( 4.0, 4.0, 0.0),
            point3( 4.0, 8.0, 0.0),
            point3( 8.0, 0.0, 0.0),
            point3( 8.0, 4.0, 0.0),
            point3( 8.0, 8.0, 0.0),
            point3(12.0, 0.0, 0.0),
            point3(12.0, 4.0, 0.0),
            point3(12.0, 8.0, 0.0),
            point3(16.0, 0.0, 0.0),
            point3(16.0, 8.0, 0.0),
            point3(16.0, 4.0, 0.0),
        ] {
            assert!(subd.vertices().contains(&p));
        }
    }

    #[test]
    fn test_subd_with_boundary_2() {
        let mesh = Mesh::from_oriented_faces(
            [
                point3( 0.0,  0.0, 0.0), // 0
                point3( 8.0,  0.0, 0.0), // 1
                point3(16.0,  0.0, 0.0), // 2
                point3( 0.0,  8.0, 0.0), // 3
                point3( 8.0,  8.0, 0.0), // 4
                point3(16.0,  8.0, 0.0), // 5
                point3( 0.0, 16.0, 0.0), // 6
                point3( 8.0, 16.0, 0.0), // 7
                point3(16.0, 16.0, 0.0), // 8
            ],
            [
                [0, 1, 4, 3],
                [1, 2, 5, 4],
                [3, 4, 7, 6],
                [4, 5, 8, 7],
            ]
        ).unwrap();

        let quad = mesh.quadrangulate(|vs| Point3::centroid(vs));
        let subd = mesh.subd(false);

        assert_eq!(quad.face_indices(), subd.face_indices());
        assert_eq!(subd.vertices().len(), 25);
        assert_eq!(subd.edge_indices().len(), 40);
        assert_eq!(subd.boundary_indices().len(), 1);
        assert_eq!(subd.boundary_indices()[0].len(), 16);

        for p in [
            point3( 0.0,  4.0, 0.0),
            point3( 0.0,  8.0, 0.0),
            point3( 0.0, 12.0, 0.0),
            point3( 1.0,  1.0, 0.0),
            point3( 1.0, 15.0, 0.0),
            point3( 4.0,  0.0, 0.0),
            point3( 4.0,  4.0, 0.0),
            point3( 4.0,  8.0, 0.0),
            point3( 4.0, 12.0, 0.0),
            point3( 4.0, 16.0, 0.0),
            point3( 8.0,  0.0, 0.0),
            point3( 8.0,  4.0, 0.0),
            point3( 8.0,  8.0, 0.0),
            point3( 8.0, 12.0, 0.0),
            point3( 8.0, 16.0, 0.0),
            point3(12.0,  0.0, 0.0),
            point3(12.0,  4.0, 0.0),
            point3(12.0,  8.0, 0.0),
            point3(12.0, 12.0, 0.0),
            point3(12.0, 16.0, 0.0),
            point3(15.0,  1.0, 0.0),
            point3(15.0, 15.0, 0.0),
            point3(16.0,  4.0, 0.0),
            point3(16.0,  8.0, 0.0),
            point3(16.0, 12.0, 0.0),
        ] {
            assert!(subd.vertices().contains(&p));
        }
    }
}
