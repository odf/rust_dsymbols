use nom::branch::alt;
use nom::character::complete::{char, digit1, space1, space0};
use nom::multi::separated_list1;
use nom::sequence::{separated_pair, tuple};
use nom::{IResult, Finish};
use nom::combinator::{map, map_opt};


pub struct DSymSpec {

}


fn integer(input: &str) -> IResult<&str, u32> {
    map_opt(digit1, map_integer)(input)
}


fn counts(input: &str) -> IResult<&str, (u32, u32)> {
    separated_pair(integer, char('.'), integer)(input)
}


fn extents(input: &str) -> IResult<&str, (u32, u32)> {
    alt((
        separated_pair(integer, space1, integer),
        map(integer, |n| (n, 2))
    ))(input)
}


fn int_list(input: &str) -> IResult<&str, Vec<u32>> {
    separated_list1(space1, integer)(input)
}


fn int_lists(input: &str) -> IResult<&str, Vec<Vec<u32>>> {
    separated_list1(tuple((space0, char(','), space0)), int_list)(input)
}


fn dsymbol(input: &str)
    -> IResult<&str, ((u32, u32), (u32, u32), Vec<Vec<u32>>, Vec<Vec<u32>>)>
{
    map(
        tuple((
            tuple((space0, char('<'), space0)),
            counts,
            tuple((space0, char(':'), space0)),
            extents,
            tuple((space0, char(':'), space0)),
            int_lists,
            tuple((space0, char(':'), space0)),
            int_lists,
            tuple((space0, char('>'), space0)),
        )),
        |(_, c, _, e, _, o, _, m, _)| (c, e, o, m)
    )(input)
}


fn map_integer(digits: &str) -> Option<u32> {
    if digits.len() <= 9 {
        Some(digits.chars().fold(0, |n, c| n * 10 + c.to_digit(10).unwrap()))
    } else {
        None
    }
}
