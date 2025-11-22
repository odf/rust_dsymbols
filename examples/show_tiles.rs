use cgmath::prelude::*;
use cgmath::{point3, vec3, vec4, Point3};
use rust_dsymbols::dsets::DSet;
use rust_dsymbols::geometry::vec_matrix::VecMatrix;
use three_d::Mat4;

use rust_dsymbols::delaney3d::pseudo_toroidal_cover;
use rust_dsymbols::display::mesh::{decompose_mesh, ItemType, Mesh};
use rust_dsymbols::dsyms::PartialDSym;
use rust_dsymbols::tilings::{Skeleton, chamber_positions, gram_matrix, invariant_basis, tile_surfaces};


fn main() {
    wrapper();
}


#[cfg(not(feature = "pprof"))]
fn wrapper() {
    run();
}


#[cfg(feature = "pprof")]
fn wrapper() {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    decompose_mesh();

    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
}


fn run() {
    // On the web, this creates a canvas instead.
    let window = three_d::Window::new(three_d::WindowSettings {
        title: "Rust 3d Test".to_string(),
        max_size: Some((1280, 720)),
        ..Default::default()
    })
    .unwrap();

    let context = window.gl();

    let mut camera = three_d::Camera::new_perspective(
        window.viewport(),
        vec3(0.0, 2.0, 8.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        three_d::degrees(25.0),
        0.1,
        1000.0,
    );

    let control = rust_dsymbols::display::controls::OrbitControl::new(1.0, 1000.0);

    //let ds_spec = "<1.1:2 3:1 2,1 2,1 2,2:3 3,4 3,4>";
    let ds_spec = "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>";

    let instances = three_d::Instances {
        transformations: vec![
            Mat4::from_scale(1.0),
        ],
        ..Default::default()
    };

    let vertex_color = three_d::Srgba::BLACK;
    let edge_color = three_d::Srgba::BLUE;
    let face_color = three_d::Srgba::RED;

    let models: Vec<_> = tiles(ds_spec).iter().flat_map(|tile_mesh|
        // todo shrink the mesh so that the individual tiles are separated
        decompose_mesh(tile_mesh).iter()
            .map(|(part_mesh, item_type)| {
                let color = match item_type {
                    ItemType::Vertex => vertex_color,
                    ItemType::Edge => edge_color,
                    ItemType::Face => face_color,
                };

                three_d::Gm::new(
                    three_d::InstancedMesh::new(
                        &context,
                        &instances,
                        &part_mesh.to_cpu_mesh()
                    ),
                    three_d::PhysicalMaterial {
                        albedo: color,
                        metallic: 0.0,
                        roughness: 0.5,
                        ..Default::default()
                    }
                )
            })
            .collect::<Vec<_>>()
        )
        .collect();

    let sun_dir = vec4(1.0, -1.0, -1.0, 0.0);

    let mut sun = three_d::DirectionalLight::new(
        &context,
        2.0,
        three_d::Srgba::WHITE,
        sun_dir.truncate()
    );

    let ambient = three_d::AmbientLight::new(
        &context,
        0.1,
        three_d::Srgba::WHITE,
    );

    window.render_loop(move |mut frame_input| {
        // This ensures a correct viewport after a window resize.
        camera.set_viewport(frame_input.viewport);

        // Camera control must be after the gui update.
        control.handle_events(&mut camera, &mut frame_input.events);

        // Moves the sun around the objects in sync with the camera.
        sun.direction = (camera.view().invert().unwrap() * sun_dir).truncate();

        frame_input.screen()
            .clear(three_d::ClearState::color_and_depth(0.8, 0.8, 0.8, 1.0, 1.0))
            .render(&camera, &models, &[&sun, &ambient]);

        // Ensures a valid return value.
        three_d::FrameOutput::default()
    });
}


fn tiles(spec: &str) -> Vec<Mesh<Point3<f64>>> {
    let ds = spec.parse::<PartialDSym>().unwrap();
    let cov = pseudo_toroidal_cover(&ds).unwrap();
    let skel = Skeleton::of(&cov);
    let basis = invariant_basis(&gram_matrix(&ds, &cov, &skel).unwrap()).transpose();
    let pos = skel.graph.vertices().iter()
        .map(|&v| (v, skel.graph.position(v)))
        .collect();

    let (scale, shift) = scale_and_shift(&cov, &skel, &basis);
    let reps = cov.orbit_reps([0, 1, 2], cov.elements());

    tile_surfaces(&cov, &skel, &pos, reps).iter().map(|(vertices, faces)| {
        let vs = vertices.iter().map(|v| {
            let v = scale * &basis * v.to_f64().unwrap() + &shift;
            point3(v[(0, 0)], v[(1, 0)], v[(2, 0)])
        });

        Mesh::from_oriented_faces(vs, faces.clone()).unwrap()
    }).collect()
}


fn scale_and_shift(cov: &PartialDSym, skel: &Skeleton, basis: &VecMatrix<f64>)
    -> (f64, VecMatrix<f64>)
{
    let ch_pos = chamber_positions(cov, skel);
    let mut s = 0.0;
    let mut n = 0;

    for d in  cov.orbit_reps([0, 2, 3], cov.elements()) {
        let v = 2.0 * basis * (&ch_pos[&(d, 1)] - &ch_pos[&(d, 0)]).to_f64().unwrap();
        s += v.norm();
        n += 1;
    };

    let scale = n as f64 / s;
    let shift = -scale * basis * ch_pos[&(1, 3)].to_f64().unwrap();

    (scale, shift)
}


trait ToCpuMesh {
    fn to_cpu_mesh(&self) -> three_d::CpuMesh;
}


impl ToCpuMesh for Mesh<Point3<f64>> {
    fn to_cpu_mesh(&self) -> three_d::CpuMesh {
        let trimesh = self.triangulate().unwrap();

        let positions: Vec<_> = trimesh.vertices().iter()
            .map(|&p| p.to_vec())
            .collect();

        let indices: Vec<_> = trimesh.face_indices().iter()
            .flatten()
            .map(|&x| x as u32)
            .collect();

        let mut mesh = three_d::CpuMesh {
            positions: three_d::Positions::F64(positions),
            indices: three_d::Indices::U32(indices),
            ..Default::default()
        };
        mesh.compute_normals();
        mesh
    }
}
