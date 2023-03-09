use crate::{common::*, constants::*, Inputs, OptionType};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
impl Inputs {
    /// Calculates the price of the option.
    /// # Requires
    /// f, k, r, q, t, sigma.
    /// # Returns
    /// f32 of the price of the option.
    /// # Example
    /// ```
    /// use black76::{Inputs, OptionType, Pricing};
    /// let inputs = Inputs::new(OptionType::Call, 100.0, 100.0, None, 0.05, 20.0/365.25, Some(0.2));
    /// let price = inputs.calc_price().unwrap();
    /// ```
    pub fn calc_price(&self) -> Result<f32, String> {
        // Calculates the price of the option
        let (nd1, nd2): (f32, f32) = calc_nd1nd2(&self)?;
        let price: f32 = match self.option_type {
            OptionType::Call => f32::max(
                0.0,
                E.powf(-self.r * self.t) * (nd1 * self.f - nd2 * self.k),
            ),
            OptionType::Put => f32::max(
                0.0,
                E.powf(-self.r * self.t) * (nd2 * self.k - nd1 * self.f),
            ),
        };
        Ok(price)
    }
}
