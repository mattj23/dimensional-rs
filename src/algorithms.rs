use ncollide2d::na::RealField;

///
pub fn preceding_index_search<N: RealField + Copy>(slice: &[N], test_value: N) -> usize {
    if slice.len() <= 1 || slice[1] > test_value {
        return 0;
    }

    let mut a = 1;
    let mut b = slice.len() - 1;
    if slice[b] <= test_value {
        return b;
    }

    while b > a + 1 {
        let check = (a + b) / 2;
        if test_value >= slice[check] {
            a = check;
        } else {
            b = check;
        }
    }
    a
}



#[cfg(test)]
mod tests {
    use super::*;
    use test_case::test_case;
    use rand::prelude::*;

    fn naive(slice: &[f64], test_value: f64) -> usize {
        if slice.len() <= 1 || slice[1] > test_value {
            return 0;
        }

        if slice[slice.len() - 1] <= test_value {
            return slice.len() - 1;
        }

        for (i, v) in slice.iter().skip(1).enumerate() {
            if *v > test_value {
                return i;
            }
        }

        slice.len() - 1
    }

    #[test_case(0, -1.0)]
    #[test_case(0, 0.05)]
    #[test_case(1, 0.1)]
    #[test_case(2, 0.25)]
    #[test_case(4, 0.5)]
    fn test_naive_search(e: usize, v: f64) {
        let test = [0.0, 0.1, 0.2, 0.3, 0.4];
        assert_eq!(e, naive(&test, v));
    }

    #[test_case(0, -1.0)]
    #[test_case(0, 0.05)]
    #[test_case(1, 0.1)]
    #[test_case(2, 0.25)]
    #[test_case(4, 0.5)]
    fn test_simple_binary_search(e: usize, v: f64) {
        let test = [0.0, 0.1, 0.2, 0.3, 0.4];
        assert_eq!(e, preceding_index_search(&test, v));
    }

    #[test]
    fn test_binary_search_random() {
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let count: usize = rng.gen_range(2..200);
            let mut values: Vec<f64> = (0..count).into_iter().map(|_| rng.gen_range(-10.0..10.0)).collect();
            values.sort_by(|a, b| a.partial_cmp(b).unwrap());

            for _ in 0..100 {
                let test = rng.gen_range(-11.0..11.0);
                let r0 = naive(&values, test);
                let r1 = preceding_index_search(&values, test);
                assert_eq!(r0, r1);
            }
        }
    }

}