use cgmath::prelude::*;
use three_d::{degrees, Camera, Event, MouseButton};

///
/// A control that makes the camera orbit around its target point.
///
#[derive(Clone, Copy, Debug)]
pub struct OrbitControl {
    /// The minimum distance to the camera target.
    pub min_distance: f32,
    /// The maximum distance to the camera target.
    pub max_distance: f32,
}

impl OrbitControl {
    pub fn new(min_distance: f32, max_distance: f32) -> Self {
        Self {
            min_distance,
            max_distance,
        }
    }

    /// Handles the events. Must be called each frame.
    pub fn handle_events(&self, camera: &mut Camera, events: &mut [Event]) -> bool {
        let target = camera.target();
        let mut change = false;
        for event in events.iter_mut() {
            match event {
                Event::MouseMotion { delta, button, handled, .. } => {
                    if let Some(button) = button {
                        if !*handled {
                            if self.apply_mouse_motion(camera, *delta, *button) {
                                *handled = true;
                                change = true;
                            }
                        }
                    }
                }
                Event::MouseWheel { delta, handled, .. } => {
                    if !*handled {
                        let dist = target.distance(camera.position());
                        self.apply_zoom(camera, delta.1, 0.001 * dist);
                        *handled = true;
                        change = true;
                    }
                }
                _ => {}
            }
        }
        change
    }

    fn apply_mouse_motion(
        &self, camera: &mut Camera, delta: (f32, f32), button: MouseButton
    ) -> bool
    {
        let (x, y) = delta;
        let target = camera.target();

        match button {
            MouseButton::Left => {
                let speed = 0.1;
                camera.rotate_around(target, speed * x, speed * y);
                true
            }
            MouseButton::Middle => {
                let speed = 0.01;
                let shift = -camera.right_direction() * x + camera.up_orthogonal() * y;
                camera.translate(speed * shift);
                true
            }
            MouseButton::Right => {
                let speed = 0.1;
                camera.roll(degrees(speed * x));
                true
            }
        }
    }

    fn apply_zoom(&self, camera: &mut Camera, delta: f32, speed: f32) {
        let target = camera.target();
        camera.zoom_towards(
            target, speed * delta, self.min_distance, self.max_distance,
        );
    }
}
